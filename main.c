/*
 * Codes written by Hang-Hyun Jo since 2025. 1. 14.
 * */

#include "head.h"

int main(int argv, char *argc[])
{
	time_t t_start, t_end, t_duration;
	long i, j, k, T, ENS0, ENS1, ens, seeds;
    long *timings, num_alter;
	double mu, nu, rho;

	char output[100], postfix[100], filename[500];
	FILE *time_out, *timings_out, *stat_out;

	ENS0 = atoi(argc[1]);
	ENS1 = atoi(argc[2]);
	mu = atof(argc[3]); // avg IET
	nu = atof(argc[4]); // avg burst size
	rho = atof(argc[5]); // correlation parameter
	T = atoi(argc[6]); // total time perid

    sprintf(postfix, "_mu%.0lf_nu%.0lf_rho%.3lf_T%ld", mu, nu, rho, T);
    printf("%s\n", postfix);

	time_out=fopen("time_record.txt", "a"); time(&t_start);
	fprintf(time_out, "\nexe_filename=%s\nstart time=%s", postfix, ctime(&t_start));
	printf("exe_filename=%s\nstart time=%s", postfix, ctime(&t_start));

    num_alter = (long)((float)T / (mu + nu) * 2.);

    timings = vector_long(0, T);

    // generate time series
    generator(timings, postfix, mu, nu, rho, T, ENS0, ENS1);

    sprintf(filename, "autocorrel%s", postfix);
    summarize_curves("./", filename, "", ENS1, 0);

    // analyze time series
    //analysis(timings, postfix, mu, nu, rho, T, num_alter, ENS0, ENS1);

	time(&t_end); t_duration=t_end-t_start;
	fprintf(time_out, "  end time=%s  duration=%ld\n", ctime(&t_end), t_duration);
	printf("  end time=%s  duration=%ld\n", ctime(&t_end), t_duration);
	fclose(time_out);

    free_vector_long(timings, 0, T);
}
