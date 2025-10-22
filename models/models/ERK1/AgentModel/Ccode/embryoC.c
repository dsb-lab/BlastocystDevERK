#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PI 4.*atan(1.)

int main(int argc, char *argv[])
{
    int i,j,k,it,ij,ji,iseed,Nmax,Ncells,*ndiv,*nn,nntot,*cfate;
    int irun,nruns,nst,nme,nmetot;
    int Fmax,iend,Ncurrent,mN2start,stdN2start,N2start;
    int N2ablate,fate2ablate,abflag;
    double *ncount,nepi1,npre1,ndp1;
    double *fepi,*fpre,*fdp,*fepi2,*fpre2,*fdp2,*ntot,*ntot2;
    double *numepi,*numpre,*numdp,*numepi2,*numpre2,*numdp2;
    double *bepi,*bpre,*bdp,*bepi2,*bpre2,*bdp2,*bcount;
    double *tvec;
    double mi,c,kf,ka,fd,kfred,ikfred;
    double tsim,mt,dt,ri,dth,tdiv,sdiv,ct,rnd,rnd1,rnd2;
    double h,ih,h0,h0h,kfm,cm,*nextdiv,rdiv;
    double alpha,Kf,nh,mh;
    double *x,*y,*z,*vx,*vy,*vz,*Fx,*Fy,*Fz,*r,*m;
    double *xi,*yi,*zi,*vxi,*vyi,*vzi;
    double *fx,*fy,*fz,*fvx,*fvy,*fvz;
    double *fxi,*fyi,*fzi,*fvxi,*fvyi,*fvzi;
    double *n,*ni,*fn,*fni;
    double nth,gth,nth2,gth2,nmax,fp,rntd,*rn;
    double n0,rnn0,p2ablate;

    char dum[60],name1[60],name2[60],name3[60];

    int ncord(int, int, int);
    int dran_ini(int, int);
    double ran1();
    void dran_gv(double *, int);

    FILE *input,*outdat,*outdat2,*outdat3,*outdat4,*outdat5;

    /**********  INPUT DATA  **********/

    strcpy(name1,argv[1]);
    sprintf(name2,".par");
    strcat(name1,name2);
    input = fopen(name1,"r");
    fscanf(input,"%i %s",&nruns,dum);
    fscanf(input,"%lg %s",&mi,dum);
    fscanf(input,"%lg %s",&ri,dum);
    fscanf(input,"%lg %s",&c,dum);
    fscanf(input,"%lg %s",&kf,dum);
    fscanf(input,"%lg %s",&kfred,dum);
    fscanf(input,"%lg %s",&ka,dum);
    fscanf(input,"%lg %s",&alpha,dum);
    fscanf(input,"%lg %s",&Kf,dum);
    fscanf(input,"%lg %s",&nh,dum);
    fscanf(input,"%lg %s",&mh,dum);
    fscanf(input,"%lg %s",&fd,dum);
    fscanf(input,"%lg %s",&nth,dum);
    fscanf(input,"%lg %s",&gth,dum);
    fscanf(input,"%lg %s",&nth2,dum);
    fscanf(input,"%lg %s",&gth2,dum);
    fscanf(input,"%i %s",&Nmax,dum);
    fscanf(input,"%lg %s",&tdiv,dum);
    fscanf(input,"%lg %s",&sdiv,dum);
    fscanf(input,"%lg %s",&tsim,dum);
    fscanf(input,"%lg %s",&mt,dum);
    fscanf(input,"%lg %s",&dt,dum);
    fscanf(input,"%i %s",&iseed,dum);
    fscanf(input,"%lg %s",&n0,dum);
    fscanf(input,"%lg %s",&rnn0,dum);
    fscanf(input,"%lg %s",&rntd,dum);
    fscanf(input,"%i %s",&mN2start,dum);
    fscanf(input,"%i %s",&stdN2start,dum);
    fscanf(input,"%i %s",&N2ablate,dum);
    fscanf(input,"%i %s",&fate2ablate,dum);
    fscanf(input,"%lg %s",&p2ablate,dum);
    fclose(input);

    strcpy(name1,argv[1]);
    sprintf(name2,".dat");
    strcat(name1,name2);
    outdat = fopen(name1,"w");

    strcpy(name1,argv[1]);
    sprintf(name2,"_ts.dat");
    strcat(name1,name2);
    outdat2 = fopen(name1,"w");

    strcpy(name1,argv[1]);
    sprintf(name2,"_nanog.dat");
    strcat(name1,name2);
    outdat3 = fopen(name1,"w");

    strcpy(name1,argv[1]);
    sprintf(name2,"_binned.dat");
    strcat(name1,name2);
    outdat4 = fopen(name1,"w");

    strcpy(name1,argv[1]);
    sprintf(name2,"_fates.dat");
    strcat(name1,name2);
    outdat5 = fopen(name1,"w");

/**********  SIMULATION CONSTANTS  **********/
    nst=(int)(tsim/dt+0.5);
    nme=(int)(mt/dt+0.5);
    nmetot = (int)(nst/nme);

    dth = dt*0.5;
    rdiv = 1.0/pow(2.0,1.0/3.0);
    ikfred = 1.0/kfred;
    Fmax = Nmax*Nmax;
    nmax = alpha;
    if (nh==2)
      nmax = alpha/(1+1/pow(2*Kf,2*mh));

    x = (double *)calloc(Nmax,sizeof(double));
    y = (double *)calloc(Nmax,sizeof(double));
    z = (double *)calloc(Nmax,sizeof(double));
    vx = (double *)calloc(Nmax,sizeof(double));
    vy = (double *)calloc(Nmax,sizeof(double));
    vz = (double *)calloc(Nmax,sizeof(double));
    Fx = (double *)calloc(Fmax,sizeof(double));
    Fy = (double *)calloc(Fmax,sizeof(double));
    Fz = (double *)calloc(Fmax,sizeof(double));
    xi = (double *)calloc(Nmax,sizeof(double));
    yi = (double *)calloc(Nmax,sizeof(double));
    zi = (double *)calloc(Nmax,sizeof(double));
    vxi = (double *)calloc(Nmax,sizeof(double));
    vyi = (double *)calloc(Nmax,sizeof(double));
    vzi = (double *)calloc(Nmax,sizeof(double));
    fx = (double *)calloc(Nmax,sizeof(double));
    fy = (double *)calloc(Nmax,sizeof(double));
    fz = (double *)calloc(Nmax,sizeof(double));
    fvx = (double *)calloc(Nmax,sizeof(double));
    fvy = (double *)calloc(Nmax,sizeof(double));
    fvz = (double *)calloc(Nmax,sizeof(double));
    fxi = (double *)calloc(Nmax,sizeof(double));
    fyi = (double *)calloc(Nmax,sizeof(double));
    fzi = (double *)calloc(Nmax,sizeof(double));
    fvxi = (double *)calloc(Nmax,sizeof(double));
    fvyi = (double *)calloc(Nmax,sizeof(double));
    fvzi = (double *)calloc(Nmax,sizeof(double));
    ndiv = (int *)calloc(Nmax,sizeof(int));
    nextdiv = (double *)calloc(Nmax,sizeof(double));
    cfate = (int *)calloc(Nmax,sizeof(int));
    r = (double *)calloc(Nmax,sizeof(double));
    m = (double *)calloc(Nmax,sizeof(double));
    n = (double *)calloc(Nmax,sizeof(double));
    ni = (double *)calloc(Nmax,sizeof(double));
    fn = (double *)calloc(Nmax,sizeof(double));
    fni = (double *)calloc(Nmax,sizeof(double));
    nn = (int *)calloc(Fmax,sizeof(int));
    rn = (double *)calloc(Nmax,sizeof(double));
    tvec = (double *)calloc(nmetot,sizeof(double));
    fpre = (double *)calloc(nmetot,sizeof(double));
    fepi = (double *)calloc(nmetot,sizeof(double));
    fdp = (double *)calloc(nmetot,sizeof(double));
    fpre2 = (double *)calloc(nmetot,sizeof(double));
    fepi2 = (double *)calloc(nmetot,sizeof(double));
    fdp2 = (double *)calloc(nmetot,sizeof(double));
    ntot = (double *)calloc(nmetot,sizeof(double));
    ntot2 = (double *)calloc(nmetot,sizeof(double));
    ncount = (double *)calloc(nmetot,sizeof(double));
    numpre = (double *)calloc(nmetot,sizeof(double));
    numepi = (double *)calloc(nmetot,sizeof(double));
    numdp = (double *)calloc(nmetot,sizeof(double));
    numpre2 = (double *)calloc(nmetot,sizeof(double));
    numepi2 = (double *)calloc(nmetot,sizeof(double));
    numdp2 = (double *)calloc(nmetot,sizeof(double));
    bpre = (double *)calloc(Nmax,sizeof(double));
    bepi = (double *)calloc(Nmax,sizeof(double));
    bdp = (double *)calloc(Nmax,sizeof(double));
    bpre2 = (double *)calloc(Nmax,sizeof(double));
    bepi2 = (double *)calloc(Nmax,sizeof(double));
    bdp2 = (double *)calloc(Nmax,sizeof(double));
    bcount = (double *)calloc(Nmax,sizeof(double));

    if (dran_ini(iseed,Nmax)==1) return 1;

/**********  SIMULATION  **********/

    for (irun=0;irun<nruns;irun++)
    {
        for (i=0;i<Nmax;i++)
        {
            n[i] = 0;
            cfate[i] = 0;
            ndiv[i] = 0;
            nextdiv[i] = 0;
            r[i] = 0;
            m[i] = 0;
            x[i] = 0;
            y[i] = 0;
            z[i] = 0;
            vx[i] = 0;
            vy[i] = 0;
            vz[i] = 0;
        }
        for (i=0;i<Fmax;i++)
        {
            Fx[i] = 0;
            Fy[i] = 0;
            Fz[i] = 0;
            nn[i] = 0;
        }
        Ncells = 1;
        n[0] = n0*(1+rnn0*(2.0*ran1(&iseed)-1.0));
        if (n[0]>nth*nmax)
            cfate[0] = 1;
        else if (n[0]<gth*nmax)
            cfate[0] = 2;
        else
          cfate[0] = 0;

        r[0] = ri;
        m[0] = mi;
        ndiv[0] = 1;
        rnd = ran1(&iseed);
        nextdiv[0] = (ndiv[0] - sdiv + rnd*2.0*sdiv)*tdiv;
        dran_gv(rn,Ncells);
        N2start = (int)(mN2start+stdN2start*rn[0]+0.5);
        abflag = 0;

        it = 0;
        while (Ncells<N2start)
        {
            for (i=1;i<Ncells;i++)
                for (j=0;j<i;j++)
                {
                    h = pow(pow(x[i]-x[j],2.0) + pow(y[i]-y[j],2.0)
                      + pow(z[i]-z[j],2.0), 0.5);
                    h0 = r[i]+r[j];
                    ih = 1.0/h;
                    h0h = h0*ih;
                    ij=i*Ncells+j;
                    if (h<ka*h0)
                    {
                        Fx[ij] = (h0h-1)*(ka*h0h-1)*ih*(x[i]-x[j]);
                        Fy[ij] = (h0h-1)*(ka*h0h-1)*ih*(y[i]-y[j]);
                        Fz[ij] = (h0h-1)*(ka*h0h-1)*ih*(z[i]-z[j]);
                    }
                    else
                    {
                        Fx[ij] = 0;
                        Fy[ij] = 0;
                        Fz[ij] = 0;
                    }
                    ji = j*Ncells+i;
                    Fx[ji] = -Fx[ij];
                    Fy[ji] = -Fy[ij];
                    Fz[ji] = -Fz[ij];
                }
            for (i=0;i<Ncells;i++)
            {
                kfm = kf/m[i];
                cm = c/m[i];
                fx[i] = vx[i];
                fy[i] = vy[i];
                fz[i] = vz[i];
                fvx[i] = -cm*vx[i];
                fvy[i] = -cm*vy[i];
                fvz[i] = -cm*vz[i];
                for (j=0;j<Ncells;j++)
                    if (j!=i)
                    {
                        ij=i*Ncells+j;
                        fvx[i] += kfm*Fx[ij];
                        fvy[i] += kfm*Fy[ij];
                        fvz[i] += kfm*Fz[ij];
                    }
                xi[i] = x[i] + dt*fx[i];
                yi[i] = y[i] + dt*fy[i];
                zi[i] = z[i] + dt*fz[i];
                vxi[i] = vx[i] + dt*fvx[i];
                vyi[i] = vy[i] + dt*fvy[i];
                vzi[i] = vz[i] + dt*fvz[i];
            }
            for (i=1;i<Ncells;i++)
                for (j=0;j<i;j++)
                {
                    h = pow(pow(xi[i]-xi[j],2.0)+pow(yi[i]-yi[j],2.0)+pow(zi[i]-zi[j],2.0),0.5);
                    h0 = r[i]+r[j];
                    ih = 1.0/h;
                    h0h = h0*ih;
                    ij=i*Ncells+j;
                    if (h<ka*h0)
                    {
                        Fx[ij] = (h0h-1)*(ka*h0h-1)*ih*(xi[i]-xi[j]);
                        Fy[ij] = (h0h-1)*(ka*h0h-1)*ih*(yi[i]-yi[j]);
                        Fz[ij] = (h0h-1)*(ka*h0h-1)*ih*(zi[i]-zi[j]);
                    }
                    else
                    {
                        Fx[ij] = 0;
                        Fy[ij] = 0;
                        Fz[ij] = 0;
                    }
                    ji=j*Ncells+i;
                    Fx[ji] = -Fx[ij];
                    Fy[ji] = -Fy[ij];
                    Fz[ji] = -Fz[ij];
                }
            for (i=0;i<Ncells;i++)
            {
                kfm = kf/m[i];
                cm = c/m[i];
                fxi[i] = vxi[i];
                fyi[i] = vyi[i];
                fzi[i] = vzi[i];
                fvxi[i] = -cm*vxi[i];
                fvyi[i] = -cm*vyi[i];
                fvzi[i] = -cm*vzi[i];
                for (j=0;j<Ncells;j++)
                    if (j!=i)
                    {
                        ij=i*Ncells+j;
                        fvxi[i] += kfm*Fx[ij];
                        fvyi[i] += kfm*Fy[ij];
                        fvzi[i] += kfm*Fz[ij];
                    }
                x[i] = x[i] + dth*(fx[i]+fxi[i]);
                y[i] = y[i] + dth*(fy[i]+fyi[i]);
                z[i] = z[i] + dth*(fz[i]+fzi[i]);
                vx[i] = vx[i] + dth*(fvx[i]+fvxi[i]);
                vy[i] = vy[i] + dth*(fvy[i]+fvyi[i]);
                vz[i] = vz[i] + dth*(fvz[i]+fvzi[i]);
            }

            ct = it*dt;
            Ncurrent = Ncells;
            for (i=0;i<Ncurrent;i++)
            {
                if (ct>nextdiv[i])
                {
                    ndiv[i]++;
                    rnd = ran1(&iseed);
                    nextdiv[i] = (ndiv[i]-sdiv + rnd*2.0*sdiv)*tdiv;
                    ndiv[Ncells] = ndiv[i];
                    rnd = ran1(&iseed);
                    nextdiv[Ncells] = (ndiv[Ncells]-sdiv + rnd*2.0*sdiv)*tdiv;
                    rnd1 = ran1(&iseed)*2.0*PI;
                    rnd2 = ran1(&iseed)*2.0*PI;
                    x[Ncells] = x[i] + r[i]*0.5*sin(rnd1)*cos(rnd2);
                    y[Ncells] = y[i] + r[i]*0.5*sin(rnd1)*sin(rnd2);
                    z[Ncells] = z[i] + r[i]*0.5*cos(rnd1);
                    x[i] = x[i] - r[i]*0.5*sin(rnd1)*cos(rnd2);
                    y[i] = y[i] - r[i]*0.5*sin(rnd1)*sin(rnd2);
                    z[i] = z[i] - r[i]*0.5*cos(rnd1);
                    vx[Ncells] = vx[i];
                    vy[Ncells] = vy[i];
                    vz[Ncells] = vz[i];
                    r[i] = rdiv*r[i];
                    r[Ncells] = r[i];
                    m[i] = m[i]*0.5;
                    m[Ncells] = m[i];
                    n[Ncells] = n[i]*(1.0+rntd*ran1(&iseed));
                    cfate[Ncells] = cfate[i];

                    n[i] = n[i]*(1.0+rntd*ran1(&iseed));

                    Ncells++;
                }
            }

            if (it%nme == 0)
            {
              ndp1 = 0;
              nepi1 = 0;
              npre1 = 0;
              for (i=0;i<Ncells;i++)
              {
                  if (cfate[i]==0)
                      ndp1 += 1;
                  if (cfate[i]==1)
                      nepi1 += 1;
                  if (cfate[i]==2)
                      npre1 += 1;
              }
              tvec[it/nme] = ct;
              fdp[it/nme] += (double)ndp1/(double)Ncells;
              fepi[it/nme] += (double)nepi1/(double)Ncells;
              fpre[it/nme] += (double)npre1/(double)Ncells;
              fdp2[it/nme] += (double)(ndp1*ndp1)/(double)(Ncells*Ncells);
              fepi2[it/nme] += (double)(nepi1*nepi1)/(double)(Ncells*Ncells);
              fpre2[it/nme] += (double)(npre1*npre1)/(double)(Ncells*Ncells);
              ntot[it/nme] += (double)(Ncells);
              ntot2[it/nme] += (double)(Ncells*Ncells);
              ncount[it/nme] += 1;
              numdp[it/nme] += (double)ndp1;
              numepi[it/nme] += (double)nepi1;
              numpre[it/nme] += (double)npre1;
              numdp2[it/nme] += (double)(ndp1*ndp1);
              numepi2[it/nme] += (double)(nepi1*nepi1);
              numpre2[it/nme] += (double)(npre1*npre1);
              bdp[Ncells] += (double)ndp1/(double)Ncells;
              bepi[Ncells] += (double)nepi1/(double)Ncells;
              bpre[Ncells] += (double)npre1/(double)Ncells;
              bdp2[Ncells] += (double)(ndp1*ndp1)/(double)(Ncells*Ncells);
              bepi2[Ncells] += (double)(nepi1*nepi1)/(double)(Ncells*Ncells);
              bpre2[Ncells] += (double)(npre1*npre1)/(double)(Ncells*Ncells);
              bcount[Ncells] += 1;
              if (irun==0)
              {
                fprintf(outdat2,"%f",ct);
                for (i=0;i<Ncells;i++)
                    fprintf(outdat2," %f",n[i]);
                fprintf(outdat2,"\n");
                fprintf(outdat5,"%f",ct);
                for (i=0;i<Ncells;i++)
                    fprintf(outdat5," %i",cfate[i]);
                fprintf(outdat5,"\n");
              }
            }

            if (Ncells == N2ablate && abflag == 0)
            {
              for (j=0;j<Ncells;j++)
                if (cfate[j]==fate2ablate || fate2ablate==3)
                  if (ran1(&iseed)<p2ablate)
                  {
                    if (nruns==1)
                      printf("Ablating cell %i with fate %i\n",j,cfate[j]);
                    cfate[j] = 3;
                  }
              for (j=0;j<Ncells;j++)
                if (cfate[j] == 3)
                {
                  for (k=j;k<Ncells-1;k++)
                  {
                    x[k] = x[k+1];
                    y[k] = y[k+1];
                    z[k] = z[k+1];
                    vx[k] = vx[k+1];
                    vy[k] = vy[k+1];
                    vz[k] = vz[k+1];
                    m[k] = m[k+1];
                    r[k] = r[k+1];
                    n[k] = n[k+1];
                    cfate[k] = cfate[k+1];
                    ndiv[k] = ndiv[k+1];
                    nextdiv[k] = nextdiv[k+1];
                  }
                  Ncells = Ncells-1;
                }
              abflag = 1;
            }

            it++;
        }

        iend = 0;
        while (iend == 0)
        {
            for (i=1;i<Ncells;i++)
                for (j=0;j<i;j++)
                {
                    h = pow(pow(x[i]-x[j],2.0) + pow(y[i]-y[j],2.0) + pow(z[i]-z[j],2.0), 0.5);
                    h0 = r[i]+r[j];
                    ih = 1.0/h;
                    h0h = h0*ih;
                    ij=i*Ncells+j;
                    if (h<ka*h0)
                    {
                        Fx[ij] = (h0h-1)*(ka*h0h-1)*ih*(x[i]-x[j]);
                        Fy[ij] = (h0h-1)*(ka*h0h-1)*ih*(y[i]-y[j]);
                        Fz[ij] = (h0h-1)*(ka*h0h-1)*ih*(z[i]-z[j]);
                    }
                    else
                    {
                        Fx[ij] = 0;
                        Fy[ij] = 0;
                        Fz[ij] = 0;
                    }
                    ji=j*Ncells+i;
                    Fx[ji] = -Fx[ij];
                    Fy[ji] = -Fy[ij];
                    Fz[ji] = -Fz[ij];
                    if (h<fd*h0)
                        nn[ij] = 1;
                    else
                        nn[ij] = 0;
                    ji=j*Ncells+i;
                    nn[ji] = nn[ij];
                }
           for (i=0;i<Ncells;i++)
            {
                kfm = kf/m[i];
                cm = c/m[i];
                fx[i] = vx[i];
                fy[i] = vy[i];
                fz[i] = vz[i];
                fvx[i] = -cm*vx[i];
                fvy[i] = -cm*vy[i];
                fvz[i] = -cm*vz[i];
                fp = 0;
                nntot = 0;
                for (j=0;j<Ncells;j++)
                    if (j!=i)
                    {
                        ij=i*Ncells+j;
                        if (cfate[i] == cfate[j])
                        {
                            fvx[i] += kfm*Fx[ij];
                            fvy[i] += kfm*Fy[ij];
                            fvz[i] += kfm*Fz[ij];
                        }
                        else
                        {
                            fvx[i] += ikfred*kfm*Fx[ij];
                            fvy[i] += ikfred*kfm*Fy[ij];
                            fvz[i] += ikfred*kfm*Fz[ij];
                        }
                        if (nn[ij]!=0)
                        {
                            fp += n[j];
                            nntot += 1;
                        }
                    }
                xi[i] = x[i] + dt*fx[i];
                yi[i] = y[i] + dt*fy[i];
                zi[i] = z[i] + dt*fz[i];
                vxi[i] = vx[i] + dt*fvx[i];
                vyi[i] = vy[i] + dt*fvy[i];
                vzi[i] = vz[i] + dt*fvz[i];

                fn[i] = alpha*pow(1+pow(n[i],nh),mh)/(pow(1+pow(n[i],nh),mh)+pow(fp/(Kf*nntot),2*mh))-n[i];
                ni[i] = n[i] + dt*fn[i];
            }

          for (i=1;i<Ncells;i++)
                for (j=0;j<i;j++)
                {
                    h = pow(pow(xi[i]-xi[j],2.0)+pow(yi[i]-yi[j],2.0)+pow(zi[i]-zi[j],2.0),0.5);
                    h0 = r[i]+r[j];
                    ih = 1.0/h;
                    h0h = h0*ih;
                    ij=i*Ncells+j;
                    if (h<ka*h0)
                    {
                        Fx[ij] = (h0h-1)*(ka*h0h-1)*ih*(xi[i]-xi[j]);
                        Fy[ij] = (h0h-1)*(ka*h0h-1)*ih*(yi[i]-yi[j]);
                        Fz[ij] = (h0h-1)*(ka*h0h-1)*ih*(zi[i]-zi[j]);
                    }
                    else
                    {
                        Fx[ij] = 0;
                        Fy[ij] = 0;
                        Fz[ij] = 0;
                    }
                    ji=j*Ncells+i;
                    Fx[ji] = -Fx[ij];
                    Fy[ji] = -Fy[ij];
                    Fz[ji] = -Fz[ij];
                    if (h<fd*h0)
                        nn[ij] = 1;
                    else
                        nn[ij] = 0;
                    ji=j*Ncells+i;
                    nn[ji] = nn[ij];
                }
           for (i=0;i<Ncells;i++)
            {
                kfm = kf/m[i];
                cm = c/m[i];
                fxi[i] = vxi[i];
                fyi[i] = vyi[i];
                fzi[i] = vzi[i];
                fvxi[i] = -cm*vxi[i];
                fvyi[i] = -cm*vyi[i];
                fvzi[i] = -cm*vzi[i];
                fp = 0;
                nntot = 0;
                for (j=0;j<Ncells;j++)
                    if (j!=i)
                    {
                        ij=i*Ncells+j;
                        if (cfate[i] == cfate[j])
                        {
                            fvxi[i] += kfm*Fx[ij];
                            fvyi[i] += kfm*Fy[ij];
                            fvzi[i] += kfm*Fz[ij];
                        }
                        else
                        {
                            fvxi[i] += ikfred*kfm*Fx[ij];
                            fvyi[i] += ikfred*kfm*Fy[ij];
                            fvzi[i] += ikfred*kfm*Fz[ij];
                        }
                        if (nn[ij]!=0)
                        {
                            fp += ni[j];
                            nntot += 1;
                        }
                    }
                x[i] = x[i] + dth*(fx[i]+fxi[i]);
                y[i] = y[i] + dth*(fy[i]+fyi[i]);
                z[i] = z[i] + dth*(fz[i]+fzi[i]);
                vx[i] = vx[i] + dth*(fvx[i]+fvxi[i]);
                vy[i] = vy[i] + dth*(fvy[i]+fvyi[i]);
                vz[i] = vz[i] + dth*(fvz[i]+fvzi[i]);

                fni[i] = alpha*pow(1+pow(ni[i],nh),mh)/(pow(1+pow(ni[i],nh),mh)+pow(fp/(Kf*nntot),2*mh))-ni[i];
                if (n[i]>gth2*nmax && n[i]<nth2*nmax)
                  n[i] = n[i] + dth*(fn[i]+fni[i]);

                if (n[i]>nth*nmax && cfate[i] == 0)
                    cfate[i] = 1;
                else if (n[i]<gth*nmax && cfate[i] == 0)
                    cfate[i] = 2;
            }

            ct = it*dt;
            Ncurrent = Ncells;
            for (i=0;i<Ncurrent;i++)
            {
                if (ct>nextdiv[i] && Ncells<Nmax-1)
                {
                    ndiv[i]++;
                    rnd = ran1(&iseed);
                    nextdiv[i] = (ndiv[i]-sdiv + rnd*2.0*sdiv)*tdiv;
                    ndiv[Ncells] = ndiv[i];
                    rnd = ran1(&iseed);
                    nextdiv[Ncells] = (ndiv[Ncells]-sdiv + rnd*2.0*sdiv)*tdiv;
                    rnd1 = ran1(&iseed)*2.0*PI;
                    rnd2 = ran1(&iseed)*2.0*PI;
                    x[Ncells] = x[i] + r[i]*0.5*sin(rnd1)*cos(rnd2);
                    y[Ncells] = y[i] + r[i]*0.5*sin(rnd1)*sin(rnd2);
                    z[Ncells] = z[i] + r[i]*0.5*cos(rnd1);
                    x[i] = x[i] - r[i]*0.5*sin(rnd1)*cos(rnd2);
                    y[i] = y[i] - r[i]*0.5*sin(rnd1)*sin(rnd2);
                    z[i] = z[i] - r[i]*0.5*cos(rnd1);
                    vx[Ncells] = vx[i];
                    vy[Ncells] = vy[i];
                    vz[Ncells] = vz[i];
                    r[i] = rdiv*r[i];
                    r[Ncells] = r[i];
                    m[i] = m[i]*0.5;
                    m[Ncells] = m[i];

                    n[Ncells] = n[i]*(1.0+rntd*ran1(&iseed));
                    cfate[Ncells] = cfate[i];

                    n[i] = n[i]*(1.0+rntd*ran1(&iseed));

                    Ncells++;
                }
            }

            if (Ncells == N2ablate && abflag == 0)
            {
              for (j=0;j<Ncells;j++)
                if (cfate[j]==fate2ablate || fate2ablate==3)
                  if (ran1(&iseed)<p2ablate)
                  {
                    if (nruns==1)
                      printf("Ablating cell %i with fate %i\n",j,cfate[j]);
                    cfate[j] = 3;
                  }
                for (j=0;j<Ncells;j++)
                  if (cfate[j] == 3)
                  {
                    for (k=j;k<Ncells-1;k++)
                    {
                      x[k] = x[k+1];
                      y[k] = y[k+1];
                      z[k] = z[k+1];
                      vx[k] = vx[k+1];
                      vy[k] = vy[k+1];
                      vz[k] = vz[k+1];
                      m[k] = m[k+1];
                      r[k] = r[k+1];
                      n[k] = n[k+1];
                      cfate[k] = cfate[k+1];
                      ndiv[k] = ndiv[k+1];
                      nextdiv[k] = nextdiv[k+1];
                    }
                    Ncells = Ncells-1;
                    j=j-1;
                  }
              abflag = 1;
            }

            if (it%nme == 0)
            {
              ndp1 = 0;
              npre1 = 0;
              nepi1 = 0;
              for (i=0;i<Ncells;i++)
              {
                if (cfate[i]==0)
                    ndp1 += 1;
                if (cfate[i]==1)
                    nepi1 += 1;
                if (cfate[i]==2)
                    npre1 += 1;
              }
              tvec[it/nme] = ct;
              fdp[it/nme] += (double)ndp1/(double)Ncells;
              fepi[it/nme] += (double)nepi1/(double)Ncells;
              fpre[it/nme] += (double)npre1/(double)Ncells;
              fdp2[it/nme] += (double)(ndp1*ndp1)/(double)(Ncells*Ncells);
              fepi2[it/nme] += (double)(nepi1*nepi1)/(double)(Ncells*Ncells);
              fpre2[it/nme] += (double)(npre1*npre1)/(double)(Ncells*Ncells);
              ntot[it/nme] += (double)(Ncells);
              ntot2[it/nme] += (double)(Ncells*Ncells);
              ncount[it/nme] += 1;
              numdp[it/nme] += (double)ndp1;
              numepi[it/nme] += (double)nepi1;
              numpre[it/nme] += (double)npre1;
              numdp2[it/nme] += (double)(ndp1*ndp1);
              numepi2[it/nme] += (double)(nepi1*nepi1);
              numpre2[it/nme] += (double)(npre1*npre1);
              bdp[Ncells] += (double)ndp1/(double)Ncells;
              bepi[Ncells] += (double)nepi1/(double)Ncells;
              bpre[Ncells] += (double)npre1/(double)Ncells;
              bdp2[Ncells] += (double)(ndp1*ndp1)/(double)(Ncells*Ncells);
              bepi2[Ncells] += (double)(nepi1*nepi1)/(double)(Ncells*Ncells);
              bpre2[Ncells] += (double)(npre1*npre1)/(double)(Ncells*Ncells);
              bcount[Ncells] += 1;
              if (irun==0)
              {
                fprintf(outdat2,"%f",ct);
                for (i=0;i<Ncells;i++)
                    fprintf(outdat2," %f",n[i]);
                fprintf(outdat2,"\n");
                fprintf(outdat5,"%f",ct);
                for (i=0;i<Ncells;i++)
                    fprintf(outdat5," %i",cfate[i]);
                fprintf(outdat5,"\n");
              }
            }

            it++;
            if (it>nst || Ncells>=Nmax-1)
            {
                iend = 1;
                bdp[Ncells] += (double)ndp1/(double)Ncells;
                bepi[Ncells] += (double)nepi1/(double)Ncells;
                bpre[Ncells] += (double)npre1/(double)Ncells;
                bdp2[Ncells] += (double)(ndp1*ndp1)/(double)(Ncells*Ncells);
                bepi2[Ncells] += (double)(nepi1*nepi1)/(double)(Ncells*Ncells);
                bpre2[Ncells] += (double)(npre1*npre1)/(double)(Ncells*Ncells);
                bcount[Ncells] += 1;
            }
        }
    }
    for (i=0;i<nmetot;i++)
        if (i>0 && tvec[i]>0 && ncount[i]>0)
        {
          fprintf(outdat,"%f %f %f %f %f %f %f %f %f\n",tvec[i],
            fdp[i]/ncount[i],fepi[i]/ncount[i],
            fpre[i]/ncount[i],fdp2[i]/ncount[i],
            fepi2[i]/ncount[i],fpre2[i]/ncount[i],ntot[i]/ncount[i],ntot2[i]/ncount[i]);
          }

    for (i=1;i<nmetot;i++)
    {
        fprintf(outdat3,"%f %f %f %f %f %f %f\n",tvec[i],
          numdp[i]/ncount[i],numepi[i]/ncount[i],
          numpre[i]/ncount[i],numdp2[i]/ncount[i],
          numepi2[i]/ncount[i],numpre2[i]/ncount[i]);
    }

    for (i=1;i<Nmax;i++)
    {
        fprintf(outdat4,"%i %f %f %f %f %f %f\n",i,
          bdp[i]/bcount[i],bepi[i]/bcount[i],
          bpre[i]/bcount[i],bdp2[i]/bcount[i],
          bepi2[i]/bcount[i],bpre2[i]/bcount[i]);
    }

    fclose(outdat);
    fclose(outdat2);
    fclose(outdat3);
    fclose(outdat4);
    fclose(outdat5);
}
