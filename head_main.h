// draw a number
double get_exponential(double xmin, double xmean){
    double rand;
    rand = 1. - genrand_real2();
    return (xmin - (xmean - xmin) * log(rand));
}

// calculate autocorrelation function
void get_autocorrel_linear(long *timings, long *timings_cumul, long T, char *postfix){
    long i, t, td, td_min, td_max, T_td, td_bin;
    double xx, n1, n2, mu1, mu2, sigma1, sigma2, autocorrel;
    char output[100], temp[100];
    FILE *autocorrel_out;

    sprintf(output, "autocorrel%s.txt", postfix); autocorrel_out = fopen(output, "w");

    timings_cumul[0] = 1; // the first event
    for(t = 1; t < T; t ++) timings_cumul[t] = timings_cumul[t - 1] + timings[t];

    td_min = 1;
    td_max = 100;
    td_bin = 1;

    fprintf(autocorrel_out, "0 1\n");
    for(td = td_min; td <= td_max; td += td_bin){
        T_td = T - td;
        n1 = timings_cumul[T_td - 1];
        n2 = timings_cumul[T - 1] - timings_cumul[td - 1];
        mu1 = n1 / (double)T_td;
        mu2 = n2 / (double)T_td;
        sigma1 = sqrt(mu1 - mu1 * mu1);
        sigma2 = sqrt(mu2 - mu2 * mu2);

        xx = 0;
        for(t = 0; t < T_td; t ++)
            if(timings[t] && timings[t + td]) xx ++;

         autocorrel = (xx / (double)T_td - mu1 * mu2) / sigma1 / sigma2;
         fprintf(autocorrel_out, "%ld %.10f\n", td, autocorrel);
     }
     fclose(autocorrel_out);
}

void get_timeseries_rho0(long *timings, long T, double mu, double nu){
    long t, i, b, iet;
    double b_real, iet_real;

    for(t = 0; t < T; t ++) timings[t] = 0;

    t = 0;
    while(1){
        b_real = get_exponential(1, nu);
        //b = (long)(b_real);
        b = (long)(b_real + 0.5);
        if(t + b < T){
            for(i = t; i < t + b; i ++) timings[i] = 1;
            t += b;
        }
        else{
            for(i = t; i < T; i ++) timings[i] = 1;
            break;
        }

        iet_real = get_exponential(1, mu);
        //iet = (long)(iet_real);
        iet = (long)(iet_real + 0.5);
        t += iet;
        if(t >= T) break;
    }
}

void get_timeseries_rho(long *timings, long T, double mu, double nu, double rho){
    long t, i, b, iet;
    double b_real, iet_real, b0_real, ci, x;

    for(t = 0; t < T; t ++) timings[t] = 0;

    b0_real = get_exponential(1, nu);
    b = (long)(b0_real + 0.5);
    for(i = 0; i < b; i ++) timings[i] = 1;
    t = b;
    while(1){
        iet_real = get_exponential(1, mu);
        iet = (long)(iet_real + 0.5);
        t += iet;
        if(t >= T) break;

        ci = rho * (1. - 2. * exp(-b0_real / nu));
        x = genrand_real2();
        b_real = nu * log(2. * ci / (ci + 1. - sqrt((ci + 1.) * (ci + 1.) - 4. * ci * x)));

        b = (long)(b_real + 0.5);
        if(t + b < T){
            for(i = t; i < t + b; i ++) timings[i] = 1;
            t += b;
        }
        else{
            for(i = t; i < T; i ++) timings[i] = 1;
            break;
        }
        b0_real = b_real;
    }
}

void generator(long *timings, char *postfix, double mu, double nu, double rho, long T, long ENS0, long ENS1){
    long ens;
    long *timings_cumul;
    char postfix1[500];

    timings_cumul = vector_long(0, T - 1);

    for(ens = ENS0; ens < ENS1; ens ++){
        printf("ens=%ld\n", ens);
	    init_genrand(10000 + ens);
	    sprintf(postfix1, "%s_ens%ld", postfix, ens); 

        if(rho == 0) get_timeseries_rho0(timings, T, mu, nu);
        else get_timeseries_rho(timings, T, mu, nu, rho);
        
        get_autocorrel_linear(timings, timings_cumul, T, postfix1);
    }

    free_vector_long(timings_cumul, 0, T - 1);
}

long count_line(char *input){
    long i;
    char c;
    FILE *file_in; 
    
    file_in = fopen(input, "r");
    i = 0;
    c = fgetc(file_in);
    while(c != EOF){
        if(c == '\n') i ++; 
        c = fgetc(file_in);
    }   
    fclose(file_in); 
    
    return i;
}

void find_moments(double **curves, long index, long ENS, double *avg, double *std){
    long i, ens, index_x;
    double x, avg0, std0;

    avg0 = 0.;
    std0 = 0.;
    for(ens = 0; ens < ENS; ens ++){
        x = curves[index][ens + 1];
        avg0 += x;
        std0 += x * x;
    }
    avg0 /= (double)ENS;
    std0 = sqrt(std0 / (double)ENS - avg0 * avg0);

    *avg = avg0;
    *std = std0;
}

// summarize curves
void summarize_curves(char *folder, char *filename, char *prefix, long ENS, long option){
    long i, j, k, x_min, x_max, ens;
    double x, y, z, **curves, avg, std;
    char input[500], output[500], postfix[100];
    FILE *file_in, *file_out;

    if(option == 0) sprintf(postfix, "");
    else if(option == 1) sprintf(postfix, "Log");

    x_min = 0;
    x_max = 0;
    for(ens = 0; ens < ENS; ens ++){
        sprintf(input, "%s%s%s_ens%ld%s.txt", folder, prefix, filename, ens, postfix);
        k = count_line(input);
        if(k > x_max) x_max = k;
    }
    printf("max line count=%ld\n", k);

    curves = matrix_double(x_min, x_max, 0, ENS);

    for(i = x_min; i <= x_max; i ++){
        for(j = 0; j <= ENS; j ++) curves[i][j] = 0;
    }

    for(ens = 0; ens < ENS; ens ++){
        sprintf(input, "%s%s%s_ens%ld%s.txt", folder, prefix, filename, ens, postfix);
        file_in = fopen(input, "r");
        k = 0;
        while(fscanf(file_in, "%lf %lf", &x, &y) && !feof(file_in)){
            curves[k][0] = x;
            curves[k][ens + 1] = y;
            k ++;
        }
        fclose(file_in);
    }

    sprintf(output, "%s%s%s.txt", folder, prefix, filename);
    file_out = fopen(output, "w");
    for(i = 0; i < k; i ++){
        find_moments(curves, i, ENS, &avg, &std);
        fprintf(file_out, "%lf %.10lf %.10lf\n", curves[i][0], avg, std);
    }
    fclose(file_out);

    free_matrix_double(curves, x_min, x_max, 0, ENS);
}
