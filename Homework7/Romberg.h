#pragma once

extern int failed_to_meet_eps_cnt;
extern int total_called_cnt;

double Romberg(double a, double b, double eps, int M, double (*f)(double x));
