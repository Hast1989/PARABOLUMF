#include <iostream>
#include <fstream>
#include<cmath>
#include <iomanip>
#include<string>
double cp, rho,sigma,alph,bett,gamma,x1,x2,k1,k2,L,Tmax,t0,u0,tau,h,eps,Q,P;
int indKap,indU, indgran, indQ, indP, n,ind;
double NormV(double* x, double* y, int n)
{
    double res;
    res = 0;
    for (int i = 0; i < n; i++)
    {
        res = (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sqrt(res);
}
double UT(double x)
{
    if (indU == 0)
    {
        return u0;
    }
    if(indU == 1)
    {
        return u0 +x*(L-x) ;
    }
    return 0;
}
double UR(double t)
{
    return u0;
}
double UL(double t)
{
    if (indgran == 0)
    {
        return u0;
    }

    return  10*sqrt(t);
}
double Pt(double t)
{
    if (indP == 0)
    {
        return 0;
    }
    if (indP == 1)
    {
        if (t <= t0)
        {
            return P;
        }
        return 0;
    }
    if (indP == 2)
    {
        if (t <= t0)
        {
            return 2*P*t;
        }
        return 0;
    }
    if (indP == 3)
    {
        if (t <= t0)
        {
            return 2 * P * (t0-t);
        }
        return 0;
    }
    if (indP == 4)
    {
        if (t <= 0.5*t0)
        {
            return 2 * P * t;
        }
        if (t <= t0)
        {
            return 2 * P * (t0 - t);;
        }
        return 0;
    }
}
double Qt(double t)
{
    if (indQ == 0)
    {
        return 0;
    }
    if (indQ == 1)
    {
        if (t <= t0)
        {
            return -Q;
        }
        return 0;
    }
    if (indQ == 2)
    {
        if (t <= t0)
        {
            return -2 * Q * t;
        }
        return 0;
    }
    if (indQ == 3)
    {
        if (t <= t0)
        {
            return -2 * Q * (t0 - t);
        }
        return 0;
    }
    if (indQ == 4)
    {
        if (t <= 0.5 * t0)
        {
            return -2 * Q * t;
        }
        if (t <= t0)
        {
            return -2 * Q * (t0 - t);;
        }
        return 0;
    }
}
double Kap(double x,double y)
{
    if (indKap == 0)
    {
        if (x <= x1)
        {
            return k1;
        }
        if (x >= x2)
        {
            return k2;
        }
        return k1 * ((x - x2) / (x1 - x2)) + k2 * ((x - x1) / (x2 - x1));
    }
    else
    {
        return 0.5 * std::pow(y, 2);
    }
    if (indKap == 2)
    {
        if (x <= x1)
        {
            return k1;
        }
        if (x >= x2)
        {
            return k2;
        }
        return k1 * ((x - x2) / (x1 - x2)) + k2 * ((x - x1) / (x2 - x1));
    }
}
double ai(double x1,double x2,double y1, double y2)
{
    return 0.5*(Kap(x1, y1) + Kap(x2, y2));
}
void Progonka(int m,int n, double* L, double* D, double* U, double* rightb, double* x)  //n - количество не граничных узлов
{
    double* alpha;
    double* betta;
    double res;
    alpha = new double[n];
    betta = new double[n];
    alpha[m] = -U[m] / D[m];
    betta[m] = rightb[m] / D[m];
    alpha[n - 1] = 0;
    betta[n - 1] = 0;
    for (int i = m+1; i < n - 1; i++)
    {
        res = L[i] * alpha[i - 1] + D[i];
        alpha[i] = -U[i] / res;
        betta[i] = (rightb[i] - L[i] * betta[i - 1]) / res;
    }
    x[n - 1] = (rightb[n - 1] - L[n - 1] * betta[n - 2]) / (L[n - 1] * alpha[n - 2] + D[n - 1]);
    for (int i = n - 2; i >= m; i--)
    {
        x[i] = alpha[i] * x[i + 1] + betta[i];
    }
    delete[] alpha;
    delete[] betta;
}
void MStep0(double tl, double tn, double* x, double* yl, double* yn, int n)
{
    double* L;
    L = new double[n + 1];
    double* res;
    res = new double[n + 1];
    double* D;
    D = new double[n + 1];
    double* U;
    U = new double[n + 1];
    double* rightb;
    rightb = new double[n + 1];
    double* a;
    a = new double[n + 1];
    double v;
    v = cp * rho * h/tau;
    for (int i = 0; i < n ; i++)
    {
        a[i + 1] = ai(x[i], x[i + 1], yl[i], yl[i + 1]);
    }
    if (sigma > eps / 100)
    {
        yn[0] = UL(tn);
        yn[n] = UR(tn);
        for (int i = 1; i < n; i++)
        {
            L[i] = -sigma * a[i] / h;
            U[i] = -sigma * a[i + 1] / h;
            D[i] = v + sigma * (a[i + 1] + a[i]) / h;
            rightb[i] = v * yl[i] + (1 - sigma) * (a[i + 1] * (yl[i + 1] - yl[i]) / h - a[i] * (yl[i] - yl[i - 1]) / h);
        }
        rightb[1] = rightb[1] + sigma * a[1] * yn[0] / h;
        rightb[n-1] = rightb[n-1] + sigma * a[n] * yn[n] / h;
        Progonka(1,n, L, D, U, rightb, yn);
    }
    else
    {
        yn[0] = UL(tn);
        yn[n] = UR(tn);
        for (int i = 1; i < n; i++)
        {
            yn[i] = (1 / v) * (a[i + 1] * (yl[i + 1] - yl[i]) / h - a[i] * (yl[i] - yl[i - 1]) / h) + yl[i];
        }
    }
    
    delete[] L;
    delete[] D;
    delete[] U;
    delete[] rightb;
    delete[] a;
    delete[] res;
}
void MStep1(double tl, double tn, double* x, double* yl, double* yn, int n)
{
    double* L;
    L = new double[n + 1];
    double* res;
    res = new double[n + 1];
    double* D;
    D = new double[n + 1];
    double* U;
    U = new double[n + 1];
    double* rightb;
    rightb = new double[n + 1];
    double* a;
    a = new double[n + 1];
    double v;
    v = cp * rho * h/tau;
    //std::cout << v << std::endl;
    for (int i = 0; i < n; i++)
    {
        a[i+1] = ai(x[i], x[i + 1],yl[i], yl[i + 1]);
    }
    if (sigma > eps / 100)
    {
        yn[0] = UL(tn);
        L[n] = -sigma * a[n] / h;
        D[n] = v / 2 + sigma * a[n] / h;
        U[n] = 0;
        rightb[n] = (v / 2) * yl[n] + sigma *  Pt(tn) + (1 - sigma) * (Pt(tl) - a[n] * (yl[n] - yl[n - 1]) / h);
        for (int i = 1; i < n; i++)
        {
            L[i] = -sigma * a[i] / h;
            U[i] = -sigma * a[i + 1] / h;
            D[i] = v + sigma * (a[i + 1] + a[i]) / h;
            rightb[i] = v * yl[i] + (1 - sigma) * (a[i + 1] * (yl[i + 1] - yl[i]) / h - a[i] * (yl[i] - yl[i - 1]) / h);
        }
        rightb[1] = rightb[1] + sigma * a[1] * yn[0] / h;
        Progonka(1,n+1, L, D, U, rightb, yn); 
        yn[0] = UL(tn);
        
    }
    else
    {
        yn[0] = UL(tn);
        yn[n] = (2 / v) * (Pt(tl) - a[n] * (yl[n] - yl[n - 1]) / h) + yl[n];
        for (int i = 1; i < n; i++)
        {
            yn[i] = (1 / v) * (a[i + 1] * (yl[i + 1] - yl[i]) / h - a[i] * (yl[i] - yl[i - 1]) / h) + yl[i];
        }
    }
    
    delete[] L;
    delete[] D;
    delete[] U;
    delete[] rightb;
    delete[] a;
    delete[] res;
}
void MStep2(double tl, double tn, double* x, double* yl, double* yn, int n)
{
    double* L;
    L = new double[n + 1];
    double* D;
    D = new double[n + 1];
    double* U;
    U = new double[n + 1];
    double* rightb;
    rightb = new double[n + 1];
    double* a;
    a = new double[n + 1];
    double v;
    v = cp * rho * h/tau;
    for (int i = 0; i < n; i++)
    {
        a[i + 1] = ai(x[i], x[i + 1], yl[i], yl[i + 1]);
    }
    if (sigma > eps / 100)
    {
        yn[n] = UR(tn);
        L[0] = 0;
        D[0] = v / 2 + sigma * a[1] / h;
        U[0] = -sigma * a[1] / h;
        rightb[0] = (v / 2) * yl[0] - sigma * Qt(tn) + (1 - sigma) * (a[1] * (yl[1] - yl[0]) / h - Qt(tl));
        for (int i = 1; i < n; i++)
        {
            L[i] = -sigma * a[i] / h;
            U[i] = -sigma * a[i + 1] / h;
            D[i] = v + sigma * (a[i + 1] + a[i]) / h;
            rightb[i] = v * yl[i] + (1 - sigma) * (a[i + 1] * (yl[i + 1] - yl[i]) / h - a[i] * (yl[i] - yl[i - 1]) / h);
        }
        rightb[n - 1] = rightb[n - 1] + sigma * a[n] * yn[n] / h;
        Progonka(0,n, L, D, U, rightb, yn);
    }
    else
    {
        yn[0] = (2 / v) * (a[1] * (yl[1] - yl[0]) / h - Qt(tl)) + yl[0];
        //std::cout << Qt(tl) << std::endl;
        yn[n] = UR(tn);
        for (int i = 1; i < n; i++)
        {
            yn[i] = (1 / v) * (a[i + 1] * (yl[i + 1] - yl[i]) / h - a[i] * (yl[i] - yl[i - 1]) / h) + yl[i];
        }
    }
    delete[] L;
    delete[] D;
    delete[] U;
    delete[] rightb;
    delete[] a;
}
void MStep3(double tl,double tn,double* x, double* yl, double* yn, int n)
{
    double* L;
    L = new double[n+1];
    double* D;
    D = new double[n+1];
    double* U;
    U = new double[n+1];
    double* rightb;
    rightb = new double[n+1];
    double* a;
    a = new double[n+1];
    double v;
    v = cp * rho*h/tau;
    for (int i = 0; i < n; i++)
    {
        a[i + 1] = ai(x[i], x[i + 1], yl[i], yl[i + 1]);
    }
    //std::cout << a[1] << std::endl;
    if (sigma > eps/100)
    {
        L[0] = 0;
        D[0] = v / 2 + sigma * a[1] / h;
        U[0] = -sigma * a[1] / h;
        rightb[0] = (v / 2) * yl[0] - sigma * Qt(tn) + (1 - sigma) * (a[1] * (yl[1] - yl[0]) / h - Qt(tl));
        L[n] = -sigma * a[n] / h;
        D[n] = v / 2 + sigma * a[n] / h;
        U[n] = 0;
        rightb[n] = (v / 2) * yl[n] + sigma * Pt(tn) + (1 - sigma) * (Pt(tl) - a[n] * (yl[n] - yl[n - 1]) / h);
        for (int i = 1; i < n; i++)
        {
            L[i] = -sigma * a[i] / h;
            U[i] = -sigma * a[i + 1] / h;
            D[i] = v + sigma * (a[i + 1] + a[i]) / h;
            rightb[i] = v * yl[i] + (1 - sigma) * (a[i + 1] * (yl[i + 1] - yl[i]) / h - a[i] * (yl[i] - yl[i - 1]) / h);
        }
        Progonka(0,n+1, L, D, U, rightb, yn);
    }
    else
    {
        yn[0] = (2/v)*(a[1]*(yl[1]-yl[0]) / h-Qt(tl))+yl[0];
        yn[n] = (2 / v) * (Pt(tl) - a[n] * (yl[n] - yl[n-1]) / h) + yl[n];
        for (int i = 1; i < n; i++)
        {
            yn[i] =(1/v)*(a[i+1]*(yl[i+1]-yl[i]) / h - a[i] * (yl[i]-yl[i-1]) / h)+yl[i];
        }
    }
    //std::cout << yl[n] << std::endl;
    delete[] L;
    delete[] D;
    delete[] U;
    delete[] rightb;
    delete[] a;
}
void test16()
{
    Q = 10.0;
    P = 10.0;
    L = 1.0;
    t0 = 0.5;
    cp = 2.0;
    rho = 0.5;
    alph = 5.0;
    bett = 0.1;
    gamma = 4.0;
    u0 = 0.04;
    k1 = 0.5;
    k2 = 1.5;
    x1 = 0.25;
    x2 = 0.5;
    tau = 0.001;
    h = 0.01;
    Tmax = 2.0;
    sigma = 0.5;
    indU = 0;
    indgran = 0;
    n = int(L / h);
    double* yl;
    double* yn;
    double* yres;
    double* x;
    double tl, tn;
    yl = new double[n + 1];
    yn = new double[n + 1];
    yres = new double[n + 1];
    x = new double[n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = i * h;
        yl[i] = UT(x[i]);
    }
    tl = 0.0;
    tn = tl + tau;
    indP = 4;
    indKap = 0;
    std::ofstream ans;
    ans.open("test16.txt");
    ans << std::setprecision(15);
    ans << h << std::endl;
    while (tl <(Tmax+tau))
    {
        MStep1(tl,tn,x,yl,yn,n);
        ans << tl<<' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << yl[i] << ' ';
            yl[i] = yn[i];
        }
        ans << std::endl;
        tl = tn;
        tn = tl + tau;

    }
    ans.close();
    delete[] yl;
    delete[] yn;
    delete[] yres;
    delete[] x;
}
void test161()
{
    Q = 10.0;
    P = 10.0;
    L = 1.0;
    t0 = 0.5;
    cp = 2.0;
    rho = 0.5;
    alph = 5.0;
    bett = 0.1;
    gamma = 4.0;
    u0 = 0.04;
    k1 = 0.1;
    k2 = 1.5;
    x1 = 0.25;
    x2 = 0.5;
    h = 0.01;
    tau = 0.0001;
    Tmax = 1.0;
    sigma = 0.5;
    indU = 1;
    indgran = 0;
    n = int(L / h);
    double* yl;
    double* yn;
    double* yres;
    double* x;
    double tl, tn;
    yl = new double[n + 1];
    yn = new double[n + 1];
    yres = new double[n + 1];
    x = new double[n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = i * h;
        yl[i] = UT(x[i]);
    }
    tl = 0.0;
    tn = tl + tau;
    indP = 0;
    indKap = 0;
    std::ofstream ans;
    ans.open("test161.txt");
    ans << std::setprecision(15);
    ans << h << std::endl;
    while (tl < (Tmax + tau))
    {
        MStep0(tl, tn, x, yl, yn, n);
        
        ans << tl << ' ';
        
        for (int i = 0; i < n + 1; i++)
        {
            ans << yl[i] << ' ';
            yl[i] = yn[i];
        }
        
        ans << std::endl;
        tl = tn;
        tn = tl + tau;

    }
    ans.close();
    delete[] yl;
    delete[] yn;
    delete[] yres;
    delete[] x;
}
void test20()
{
    Q = 10.0;
    P = 10.0;
    L = 4.0;
    t0 = 0.5;
    cp = 1.0;
    rho = 0.75;
    alph = 2.0;
    bett = 1.5;
    gamma = 2.0;
    u0 = 0.5;
    k1 = 1;
    k2 = 1;
    x1 = 1;
    x2 = 2;
    h = 0.01;
    tau = 0.1;
    Tmax = 2.0;
    sigma = 0.5;
    indU = 0;
    h = 0.2;
    indgran = 0;
    n = int(L / h);
    double* yl;
    double* yn;
    double* yres;
    double* x;
    double tl, tn;
    yl = new double[n + 1];
    yn = new double[n + 1];
    yres = new double[n + 1];
    x = new double[n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = i * h;
        yl[i] = UT(x[i]);
    }
    tl = 0.0;
    tn = tl + tau;
    indQ = 4;
    indKap = 0;

    std::ofstream ans;
    ans.open("test20.txt");
    ans << std::setprecision(15);
    ans << h << std::endl;
    while (tl < (Tmax + tau))
    {
        MStep2(tl, tn, x, yl, yn, n);
        ans << tl << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << yl[i] << ' ';
            yl[i] = yn[i];
        }
        ans << std::endl;
        tl = tn;
        tn = tl + tau;

    }
    ans.close();
    delete[] yl;
    delete[] yn;
    delete[] yres;
    delete[] x;
}
void test201()
{
    Q = 10.0;
    P = 10.0;
    L = 1.0;
    t0 = 0.5;
    cp = 1.0;
    rho = 0.75;
    alph = 2.0;
    bett = 1.5;
    gamma = 2.0;
    u0 = 0.5;
    k1 = 0.5;
    k2 = 1.5;
    x1 = 0.5;
    x2 = 0.6;
    h = 0.01;
    tau = 0.0001;
    Tmax = 1.0;
    sigma = 0.5;
    indU = 1;
    indgran = 0;
    n = int(L / h);
    double* yl;
    double* yn;
    double* yres;
    double* x;
    double tl, tn;
    yl = new double[n + 1];
    yn = new double[n + 1];
    yres = new double[n + 1];
    x = new double[n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = i * h;
        yl[i] = UT(x[i]);
    }
    tl = 0.0;
    tn = tl + tau;
    indQ = 1;
    indKap = 0;
    std::ofstream ans;
    ans.open("test201.txt");
    ans << std::setprecision(15);
    ans << h << std::endl;
    while (tl < (Tmax + tau))
    {
        MStep0(tl, tn, x, yl, yn, n);
        ans << tl << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << yl[i] << ' ';
            yl[i] = yn[i];
        }
        ans << std::endl;
        tl = tn;
        tn = tl + tau;

    }
    ans.close();
    delete[] yl;
    delete[] yn;
    delete[] yres;
    delete[] x;
}
void test162()
{
    Q = 10.0;
    P = 10.0;
    L = 1.0;
    t0 = 0.5;
    cp = 2.0;
    rho = 0.5;
    alph = 5.0;
    bett = 0.1;
    gamma = 4.0;
    u0 = 0.04;
    k1 = 0.1;
    k2 = 1.5;
    x1 = 0.25;
    x2 = 0.5;
    h = 0.01;
    tau = 0.0001;
    Tmax = 2.0;
    sigma = 0.5;
    indU = 1;
    indgran = 0;
    n = int(L / h);
    double* yl;
    double* yn;
    double* yres;
    double* x;
    double tl, tn;
    yl = new double[n + 1];
    yn = new double[n + 1];
    yres = new double[n + 1];
    x = new double[n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = i * h;
        yl[i] = UT(x[i]);
    }
    tl = 0.0;
    tn = tl + tau;
    indP = 0;
    indQ = 0;
    indKap = 0;
    std::ofstream ans;
    ans.open("test162.txt");
    ans << std::setprecision(15);
    ans << h << std::endl;
    while (tl < (Tmax + tau))
    {
        MStep3(tl, tn, x, yl, yn, n);
        ans << tl << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << yl[i] << ' ';
            yl[i] = yn[i];
        }
        ans << std::endl;
        tl = tn;
        tn = tl + tau;

    }
    ans.close();
    delete[] yl;
    delete[] yn;
    delete[] yres;
    delete[] x;
}
void test202()
{
    Q = 10.0;
    P = 10.0;
    L = 1.0;
    t0 = 0.5;
    cp = 1.0;
    rho = 0.75;
    alph = 2.0;
    bett = 1.5;
    gamma = 2.0;
    u0 = 0.5;
    k1 = 0.5;
    k2 = 1.5;
    x1 = 0.5;
    x2 = 0.6;
    h = 0.01;
    tau = 0.0001;
    Tmax = 1.0;
    sigma = 0.5;
    indU = 1;
    indgran = 0;
    n = int(L / h);
    double* yl;
    double* yn;
    double* yres;
    double* x;
    double tl, tn;
    yl = new double[n + 1];
    yn = new double[n + 1];
    yres = new double[n + 1];
    x = new double[n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = i * h;
        yl[i] = UT(x[i]);
    }
    tl = 0.0;
    tn = tl + tau;
    indP = 0;
    indQ = 0;
    indKap = 0;
    std::ofstream ans;
    ans.open("test202.txt");
    ans << std::setprecision(15);
    ans << h << std::endl;
    while (tl < (Tmax + tau))
    {
        MStep3(tl, tn, x, yl, yn, n);
        ans << tl << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << yl[i] << ' ';
            yl[i] = yn[i];
        }
        ans << std::endl;
        tl = tn;
        tn = tl + tau;

    }
    ans.close();
    delete[] yl;
    delete[] yn;
    delete[] yres;
    delete[] x;
}
void test163()
{
    Q = 10.0;
    P = 10.0;
    L = 9.8;
    t0 = 0.5;
    cp = 2.0;
    rho = 0.5;
    alph = 5.0;
    bett = 0.1;
    gamma = 4.0;
    u0 = 0;
    k1 = 0.5;
    k2 = 1.5;
    x1 = 0.25;
    x2 = 0.5;
    h = 0.2;
    tau = 0.0001;
    Tmax = 2.0;
    sigma = 1;
    indU = 0;
    indgran = 1;
    n = int(L / h);
    double* yl;
    double* yn;
    double* yres1;
    double* yres2;
    double* x;
    double tl, tn;
    yl = new double[n + 1];
    yn = new double[n + 1];
    yres1 = new double[n + 1];
    yres2 = new double[n + 1];
    x = new double[n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = i * h;
        yl[i] = UT(x[i]);
    }
    tl = 0.0;
    tn = tl + tau;
    indP = 0;
    indKap = 1;
    std::ofstream ans;
    ans.open("test163.txt");
    ans << std::setprecision(15);
    ans << h << std::endl;
    
    while (tl < (Tmax + tau))
    {
        ans << tl << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << yl[i] << ' ';
        }
        ans << std::endl;
        ind = 1;
        MStep1(tl, tn, x, yl, yres1, n);
        MStep1(tl, tn, x, yres1, yres2, n);
        while (NormV(yres1, yres2, n + 1) > eps)
        {
            ind++;
            for (int i = 0; i < n + 1; i++)
            {
                yres1[i] = yres2[i];
            }
            MStep1(tl, tn, x, yres1, yres2, n);
            if (ind > 3)
            {
                break;
            }
        }
        //std::cout << ind << std::endl;
        for (int i = 0; i < n + 1; i++)
        {
            yl[i] = yres2[i];
        }
        tl = tn;
        tn = tl + tau;

    }
    ans.close();
    delete[] yl;
    delete[] yn;
    delete[] yres1;
    delete[] yres2;
    delete[] x;
}
void test203()
{
    Q = 0.01;
    P = 0.01;
    L = 4;
    t0 = 2;
    cp = 1;
    rho = 0.01;
    alph = 5.0;
    bett = 0.1;
    gamma = 4.0;
    u0 = 0;
    k1 = 1000;
    k2 = 1;
    x1 = 0.5;
    x2 = 2;
    h = 0.2;
    tau = 0.0001;
    Tmax = 5.0;
    sigma = 0.5;
    indU = 0;
    indgran = 1;
    h = 0.4;
    n = int(L / h);
    double* yl;
    double* yn;
    double* yres1;
    double* yres2;
    double* x;
    double tl, tn;
    yl = new double[n + 1];
    yn = new double[n + 1];
    yres1 = new double[n + 1];
    yres2 = new double[n + 1];
    x = new double[n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = i * h;
        yl[i] = UT(x[i]);
    }
    tl = 0.0;
    tn = tl + tau;
    indP = 4;
    indQ = 4;
    indKap = 2;
    std::ofstream ans;
    ans.open("test203е.txt");
    ans << std::setprecision(15);
    ans << h << std::endl;

        while (tl < (Tmax + tau))
        {
            MStep3(tl, tn, x, yl, yn, n);
            ans << tl << ' ';
            for (int i = 0; i < n + 1; i++)
            {
                ans << yl[i] << ' ';
                yl[i] = yn[i];
            }
            ans << std::endl;
            tl = tn;
            tn = tl + tau;

        }

    ans.close();
    delete[] yl;
    delete[] yn;
    delete[] yres1;
    delete[] yres2;
    delete[] x;
}
int main()
{
    eps = 0.0001;
    test16();
    //test203();
    //test161();
    //test162();
    //test201();
    //test202();
    //test163();
    //std::cout << n << std::endl;
    std::cout << "Hello World!\n";
}

