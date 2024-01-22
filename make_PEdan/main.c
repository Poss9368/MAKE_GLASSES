//
//  main.c
//  vidrios_2D
//
//  Created by cjvillarroel on 5/6/19.
//  Copyright © 2019 cjvillarroel. All rights reserved.
//

#include <stdio.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>

#define    N                             1024*1
#define    DOUBLE_N                      ((double)N)
#define    INV_N                         1/DOUBLE_N


#define    UNI                           ((double)rand()/((double)RAND_MAX + 1.0))
//#define    RHO                          0.875
//#define    L                            sqrt(((1.4*1.4+1)*M_PI*DOUBLE_N)/(2*RHO))
//#define    HALF_L                       L/2

double   RHO;     //                       0.9
double   L ;      //                       sqrt(((1.4*1.4+1)*M_PII*DOUBLE_N)/(2*RHO))
double   HALF_L; //                       L/2

#define    M_PII                            3.14159265358979323846264338327


#define    MAX_CELLS                        (N)


#define    MAX_NEBZ                         256
#define    DR                               0.2
#define    C_cutoff                         1.38541802488147397906623679177079
#define    C_e                              1.0


double rx[N], ry[N], fx[N], fy[N], px[N], py[N];
int size[N];
int serial=0;

int nebz[N][MAX_NEBZ]; //maximum we give here MAX_NEBZ nebz.
int numOfNebz[N]; //number of elements in nebz[N][MAX_NEBZ];
double listCutOffSqrd=(1.4*C_cutoff+DR)*(1.4*C_cutoff+DR);
double cellCutOff=(1.4*C_cutoff+DR);


double M[N*2*MAX_NEBZ][3];
double B[2*N];

int N_con_real;
int numInCell[MAX_CELLS];
int cells[MAX_CELLS][MAX_NEBZ];
int nebListCounter;

const double SizeSqrd[2][2] = {{C_cutoff*C_cutoff,1.2*C_cutoff*1.2*C_cutoff},{1.2*C_cutoff*1.2*C_cutoff,1.4*C_cutoff*1.4*C_cutoff}}; // [0][0] = small-small; [0][1] = [1][0] = large-small; [1][1] = large-large;
const double Size_e[2][2] = {{1.00*C_cutoff,1.2*C_cutoff},{1.2*C_cutoff,1.4*C_cutoff}};

double typicalForce,typicalGrad,force_neta,Vcontrol;
int em=0;
double kinetic,T;

double max_D,tota_D=0;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradiente conjugado

void impresor(void);
void def_porte(void);
void def_in_pos(void);
void updateNebzLists(void);
void updateNebzLists_malo(void);
double desplaza1(double dt);
double fireAdvanceTime(double Factor);
void fixDrift(void);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void def_porte(void)
{
    
    for (int i =0 ; i<N; i++)
    {
        if (i<DOUBLE_N/2.0)
        {
            size[i]=0;
        }
        else
        {
            size[i]=1;
        }
    }
    
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void def_in_pos(void)
{
    for (int i =0 ; i<N; i++)
    {
        rx[i]=UNI*L;
        ry[i]=UNI*L;
    }
    
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void updateNebzLists(void)
{
    int i,j,x,y,current,m,m2;
    int a,b,k,l,numHere,numThere,target,w;
    double dx,dy,r2,invCellSize;
    
    nebListCounter = 0;
    m =(int)(L/cellCutOff);    //the length of cells and numInCell...
    m2 = m*m;
    invCellSize = m/L;     // for reduced coordinates...

    
    //set number in each cell to zero
    for (i=0; i<m2; i++)
        numInCell[i] = 0;
    
    //********cells sorting********
    for (i=0;i<N;i++){
        current =  m*(int)(ry[i]*invCellSize) + (int)(rx[i]*invCellSize);
         cells[current][numInCell[current]] = i;
        numInCell[current]++;
    }
    
    //set all neb lists to zero
    for (i=0; i<N; i++)
        numOfNebz[i] = 0;
    
    //z is layer...
    for (y=0;y<m;y++){
        for (x=0;x<m;x++){
            
            current = m*y + x; //this cell index.
            
            numHere = numInCell[current];
            //first check interactions within a cell
            for (k=0;k<numHere-1;k++){
                for (l=k+1;l<numHere;l++){
                    i = cells[current][k];
                    j = cells[current][l];
                    dx = rx[j] - rx[i];
                    dy = ry[j] - ry[i];
                    r2 =( dx*dx + dy*dy);
                    if (r2 < listCutOffSqrd){
                        
                        nebz[i][numOfNebz[i]] = j;
                        nebz[j][numOfNebz[j]] = i;
                        numOfNebz[i]++;
                        numOfNebz[j]++;
                    }
                }
            }
            
            int shift[4][2] = {{0,1},{1,0},{1,1},{1,-1}};
            //now there are 4 nebz that need to be checked.
            for (w=0;w<4;w++){
                a = x + shift[w][0];
                if (a<0) a = m-1;
                else if (a==m) a = 0;
                b = y + shift[w][1];
                if (b<0) b = m-1;
                else if (b==m) b = 0;;
                
                target = a + m*b;
                numThere = numInCell[target];
                for (k=0;k<numHere;k++){
                    for (l=0;l<numThere;l++){
                        i = cells[current][k];
                        j = cells[target][l];
                        dx = rx[j] - rx[i];
                        dy = ry[j] - ry[i];
                        // mess due to periodic boundary conditions
                        if ( dx >= HALF_L )
                            dx -= L;
                        else if ( dx < -HALF_L )
                            dx += L;
                        if ( dy >= HALF_L )
                            dy -= L;
                        else if ( dy < -HALF_L )
                            dy += L;
                        // end of mess
                        r2 =( dx*dx + dy*dy);
                        if (r2 < listCutOffSqrd){
                            nebz[i][numOfNebz[i]] = j;
                            nebz[j][numOfNebz[j]] = i;
                            numOfNebz[i]++;
                            numOfNebz[j]++;
                        }
                    }
                }
            }
        }
    }
    
    return;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void updateNebzLists_malo(void)
{
    
    
    double dx,dy,r2;
    for (int i=0; i<N; i++)
        numOfNebz[i] = 0;
    
    
    for (int i =0; i<N; i++)
    {
        for (int j =0; j<N; j++)
        {
            
            if (i<j) {
                dx = rx[j] - rx[i];
                dy = ry[j] - ry[i];
                // mess due to periodic boundary conditions
                if ( dx >= HALF_L )
                    dx -= L;
                else if ( dx < -HALF_L )
                    dx += L;
                if ( dy >= HALF_L )
                    dy -= L;
                else if ( dy < -HALF_L )
                    dy += L;
                // end of mess
                r2 = ( dx*dx + dy*dy);
                if (r2 < listCutOffSqrd)
                {
                    nebz[i][numOfNebz[i]] = j;
                    nebz[j][numOfNebz[j]] = i;
                    numOfNebz[i]++;
                    numOfNebz[j]++;
                }
                
            }
        }
    }
    
    
    
    //actualiza++;
    
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double fuerza_contactos; int potencial=0;
void calculateForces(void)
{
    int i,j,m,k;
    double goo,dx,dy,r2,r,sigma;
    double temp,sigma2;
    typicalForce=0;
    
    N_con_real=0;
    fuerza_contactos=0;
    
    for (i=0;i<N;i++){
        fx[i] = 0.0; fy[i] = 0.0;
    }
    
    for (i=0;i<N;i++){
        
        m = numOfNebz[i];
        
        for (k=0;k<m;k++){
            j = nebz[i][k];
            if (j > i){
                dx = rx[j] - rx[i];
                dy = ry[j] - ry[i];
                // mess due to periodic boundary conditions
                if ( dx >= HALF_L )
                    dx -= L;
                else if ( dx < -HALF_L )
                    dx += L;
                if ( dy >= HALF_L )
                    dy -= L;
                else if ( dy < -HALF_L )
                    dy += L;
                // end of mess
                r2 = ( dx*dx + dy*dy );
                sigma2=SizeSqrd[size[i]][size[j]];
                if ( r2 < sigma2){
                    r=sqrt(r2);
                    sigma=Size_e[size[i]][size[j]];

                    if (potencial==0)
                    {
                        //goo=C_e*(1-r/sigma)*(1/sigma);   // armonico
                        goo=C_e*pow((1-r/sigma),1.5)*(1/sigma); // hertzian
                    }
                    else
                    {
                        goo=C_e*(10*sigma2*sigma2*sigma2*sigma2*sigma2/(r2*r2*r2*r2*r2*r)+0.625201*r2*r/(sigma2*sigma2)-1.4*r/sigma2);
                        
                    }
                    
                    
                    typicalForce += goo*goo;
                    
                    temp = dx*goo/r;
                    fx[j] += temp;
                    fx[i] -= temp;
                    temp = dy*goo/r;
                    fy[j] += temp;
                    fy[i] -= temp;
                    
                    N_con_real++;
                    fuerza_contactos+=goo*r;
                    
                    
                }
            }
        }
        
        
    }
    
    
    typicalForce = sqrt(typicalForce/DOUBLE_N);
    typicalGrad = 0.0;
    force_neta=0;
    double temmm=0;
    for (i=0;i<N;i++)
    {
        temmm=fx[i]*fx[i] + fy[i]*fy[i];
        typicalGrad += temmm;
        force_neta += sqrt(temmm);
        
    }
    typicalGrad = sqrt(typicalGrad/DOUBLE_N);
    Vcontrol=typicalGrad/typicalForce;
    force_neta=force_neta/DOUBLE_N;
    fuerza_contactos=fuerza_contactos/(L*L);
    
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double desplaza1(double dt)
{
    double tem;
    max_D=0;
   // double dt=0.9;
    for (int i=0; i<N; i++)
    {
        rx[i]=rx[i]+fx[i]*dt;
        if (rx[i] >= L){
            rx[i] = rx[i] - L;
        }
        else if (rx[i] < 0.0){
            rx[i] = rx[i] + L;
        }
        
        ry[i]=ry[i]+fy[i]*dt;
        if (ry[i] >= L){
            ry[i] = ry[i] - L;
        }
        else if (ry[i] < 0.0){
            ry[i] = ry[i] + L;
        }
        
        tem=fx[i]*fx[i]+fy[i]*fy[i];
        if (max_D<tem)
        {
            max_D=tem;
        }
    }
    
    max_D=sqrt(max_D)*dt;
    return max_D;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void evol_tem1(double dt,int step)
{
    for (int i =0; i<step; i++)
    {
        max_D=desplaza1(dt);
        tota_D=tota_D+max_D;
        if (tota_D>DR)
        {
            updateNebzLists();
            tota_D=0;
        }
        calculateForces();
        printf("\n %d:  contro %0.12lf force neta %0.12lf   %d ", i+1,Vcontrol,force_neta, N_con_real);
    }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define    ALPHA_START                    0.1
#define    TIME_STEP                    0.05
#define    DELTA_T_MAX                  (10.0*TIME_STEP)
#define    N_MIN                        5
#define    F_INC                        1.1
#define    F_DEC                        0.5
#define    F_ALPHA                      0.99

#define    MAX_ITERATIONS                       100000
#define    RATIO_TOL                            1e-14


double timeStep,alpha;
int numOfStepsSinceNegativePower,itr;
double fireAdvanceTime(double Factor)
{
    int i;
    double temp,dx,dy;
    double power,mod_force,inv_mod_force,mod_velocity;
    max_D = 0;
    
    for (i=0; i<N; i++){
        px[i] += 0.5*timeStep*fx[i];
        py[i] += 0.5*timeStep*fy[i];
        
        dy = timeStep*py[i]; //real space displacement
        temp = ry[i]  + dy;
        if (temp >= L){
            ry[i] = temp - L;
        }
        else if (temp < 0.0){
            ry[i] = temp + L;
            
        }
        else
            ry[i] = temp;
        
        dx = timeStep*px[i]; //real space displacement
        temp = rx[i] + dx;
        if (temp >= L){
            rx[i] = temp - L;
        }
        else if (temp < 0.0){
            rx[i] = temp + L;
            //  ghost_rx[i] += LENGTH;
        }
        else
            rx[i] = temp;
        
        temp = dx*dx + dy*dy;
        if (temp > max_D)
            max_D = temp;
    }
    max_D= sqrt(max_D);
    
    
    calculateForces();
    for (i=0; i<N; i++){
        px[i] += 0.5*timeStep*fx[i];
        py[i] += 0.5*timeStep*fy[i];
    }
    
    
    //calculate power
    power = 0.0;
    mod_force = 0.0;
    mod_velocity = 0.0;
    for (i=0;i<N;i++){
        power += fx[i]*px[i] + fy[i]*py[i] ;
        mod_force += fx[i]*fx[i] + fy[i]*fy[i];
        mod_velocity += px[i]*px[i] + py[i]*py[i];
    }
    inv_mod_force = 1.0/sqrt(mod_force);
    mod_velocity = sqrt(mod_velocity);
    
    for (i=0;i<N;i++){
        px[i] = (1.0 - alpha)*px[i] + alpha*mod_velocity*fx[i]*inv_mod_force;
        py[i] = (1.0 - alpha)*py[i] + alpha*mod_velocity*fy[i]*inv_mod_force;
    }
    
    if ( power > 0.0 ){
        if ( numOfStepsSinceNegativePower > N_MIN ){
            if ( F_INC*timeStep < DELTA_T_MAX*Factor )
                timeStep *= F_INC;
            else
                timeStep = DELTA_T_MAX*Factor;
            alpha *= F_ALPHA;
        }
        numOfStepsSinceNegativePower++;
    }
    else{
        numOfStepsSinceNegativePower = 0;
        timeStep *= F_DEC;
        for (i=0;i<N;i++){ //freeze the system
            px[i] = 0.0;
            py[i] = 0.0;
        }
        alpha = ALPHA_START;
    }
    return max_D;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int fire(double Factor)
{
    numOfStepsSinceNegativePower = 0;
    alpha = ALPHA_START;
    timeStep = TIME_STEP*Factor;
    itr = 0;
    while ( itr < MAX_ITERATIONS && force_neta > RATIO_TOL )
    {
        max_D=fireAdvanceTime(Factor);
        tota_D=tota_D+max_D;
        if (tota_D>DR)
        {
            updateNebzLists();
            tota_D=0;
        }
        
         printf("\n %d:  contro %0.12lf force neta %0.12lf   %d ", itr+1,Vcontrol,force_neta, N_con_real);
        itr++;
        
    }
    if ( itr == MAX_ITERATIONS ){
        printf("\n        failed to minimize\n\n");
        return 0;
    }
    else
    {
        return 1;
    }
    
}  // metodo de edan, el segundo en usarse.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define    TIME_STEP_TER                     0.001
#define    HALF_TIME_STEP_TER                0.5*TIME_STEP_TER
void fixDrift(){
    int i;
    double Px,Py;
    Px = 0.0; Py = 0.0;
    for (i=0;i<N;i++){
        Px += px[i];
        Py += py[i];
    }
    Px = Px/DOUBLE_N; Py = Py/DOUBLE_N;;
    for (i=0;i<N;i++){
        px[i] -= Px;
        py[i] -= Py;
    }
    return;
}
void advanceTime(double thermostat_timescale, double Factor){
    int i;
    double temp,dx,dy,inst_T=0,xi_T;
    max_D = 0;
    
    for (i=0; i<N; i++){
        px[i] += HALF_TIME_STEP_TER*fx[i]*Factor;
        py[i] += HALF_TIME_STEP_TER*fy[i]*Factor;
    
        
        dx = TIME_STEP_TER*px[i]*Factor; //real space displacement
        temp = rx[i] + dx;
        if (temp >= L)
            rx[i] = temp - L;
        else if (temp < 0.0)
            rx[i] = temp + L;
        else
            rx[i] = temp;
        
        dy = TIME_STEP_TER*py[i]*Factor; //real space displacement
        temp = ry[i] + dy;
        if (temp >= L)
            ry[i] = temp - L;
        else if (temp < 0.0)
            ry[i] = temp + L;
        else
            ry[i] = temp;
        
        temp = dx*dx + dy*dy;
        if (temp > max_D)
            max_D = temp;
    }
    
    max_D= sqrt(max_D);
    tota_D=tota_D+max_D;
    if (tota_D>DR)
    {
        updateNebzLists();
        tota_D=0;
    }
    
    calculateForces();
    kinetic = 0.0;
    for (i=0; i<N; i++){
        px[i] += HALF_TIME_STEP_TER*fx[i]*Factor;
        py[i] += HALF_TIME_STEP_TER*fy[i]*Factor;
        kinetic += px[i]*px[i] + py[i]*py[i];
    }
    kinetic = 0.5*kinetic;
    
    /************************** thermostating part *************************************/
    if ( thermostat_timescale > 0.0 )
    {
                inst_T = 2.0*kinetic/(2.0*DOUBLE_N);
                xi_T = sqrt(1.0 + TIME_STEP_TER*Factor*(T/inst_T - 1.0)/thermostat_timescale);
                for (i=0;i<N;i++)
                {
                        px[i] *= xi_T;
                        py[i] *= xi_T;
                }
    }
    /********************************************************************************/
    printf("\n %lf %0.16lf %d", inst_T, force_neta,N_con_real);
    return;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void termostado(double duration, double thermostat_timescale,double Factor)
{
    int steps,k;

    steps = (int)(duration/TIME_STEP);
    for (k=0; k<steps; k++){
            advanceTime(thermostat_timescale,Factor);
            if ( !(k%32) )
            {
                    fixDrift();
            }
    }
    return;
    
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void impresor()
{
    FILE *fp;
    char archi[150];
    
    sprintf(archi,"result/vidrio2D_HERT_%d_%lf_%05d.dat",N,RHO,serial);
    fp = fopen(archi,"w");
    
    
    fprintf(fp,"%.18f\t%.18f\t%d\n",L,RHO,N);
    
    for (int i = 0; i<N; i++)
    {
        fprintf(fp,"%.18f\t%.18f\t%d\n",rx[i],ry[i],size[i]);
        
    }
    fclose(fp);
    
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int argc, const char * argv[])
{
     int ini; int fin;
     //sscanf(argv[1],"%d",&ini);   //rho
     //sscanf(argv[2],"%d",&fin);   //rho
     //sscanf(argv[3],"%lf",&RHO);   //rho
     //L =sqrt(((1.4*1.4+1)*M_PII*DOUBLE_N)/(2*RHO));
     //HALF_L=L/2;
    
    ini=0; fin=1;
    RHO=0.975;
    //L =sqrt(((1.4*1.4+1)*M_PII*DOUBLE_N)/(2*RHO));
    L=sqrt(DOUBLE_N/RHO);
    HALF_L=L/2;
    
    for (int tt=ini; tt<fin+ini; tt++)
    {
        potencial=0;  serial=tt;
        srand( serial + 126 );
        //double ll=L;  printf("%lf   %d    largo caja: %g    area: %lf\n",RHO,serial,ll,ll*ll);
        
        
        def_porte();
        def_in_pos();
        updateNebzLists(); tota_D=0;
        calculateForces();
        //evol_tem1(0.9);
        fire(1.0);
        
        for (int i=0;i<N;i++)
        {
            px[i] =0;
            py[i] =0;
        }
        
        potencial=1;
        updateNebzLists();  tota_D=0;
        calculateForces();
        T=0; termostado(1000, 0.1,1);
        evol_tem1(0.00001,1000);

        
        /*for (int i=0;i<N;i++)
        {
            px[i] =0;
            py[i] =0;
        }
      
        updateNebzLists();
        calculateForces();
        

       //evol_tem1(0.000001);
        for (int i=0;i<N;i++)
        {
            px[i] =0;
            py[i] =0;
        }
        fire(0.0001);*/
        
        impresor();
    }

}


