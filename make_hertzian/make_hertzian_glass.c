//
//  main.c
//  vidrios_2D
//
//  Created by cjvillarroel 
//  Copyright © 2019 cjvillarroel. All rights reserved.
//

#include <stdio.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>


#define    N                            1024  // number of particles, perfect square required
#define    DOUBLE_N                     ((double)N)
#define    INV_N                        1/DOUBLE_N


#define    UNI                          ((double)rand()/((double)RAND_MAX + 1.0))

double   RHO;     //                      0.9
double   L ;      //                      sqrt(((1.4*1.4+1)*M_PII*DOUBLE_N)/(2*RHO))
double   HALF_L;  //                      L/2

#define    M_PII                         3.14159265358979323846264338327


#define    MAX_CELLS                       (N)


#define    MAX_NEBZ                         256
#define    DR                               0.2
#define    C_e                              1.0


double rx[N], ry[N], fx[N], fy[N], px[N], py[N];
int size[N];
unsigned int serial=0;

int nebz[N][MAX_NEBZ]; //maximum we give here MAX_NEBZ nebz.
int numOfNebz[N]; //number of elements in nebz[N][MAX_NEBZ];
double listCutOffSqrd=(1.4*2+DR)*(1.4*2+DR);
double cellCutOff=(1.4*2+DR);

double M[N*2*MAX_NEBZ][3];
double B[2*N];

int numInCell[MAX_CELLS];
int cells[MAX_CELLS][MAX_NEBZ];
int nebListCounter;

const double SizeSqrd[2][2] = {{4.00,2.4*2.4},{2.4*2.4,2.8*2.8}}; // [0][0] = small-small; [0][1] = [1][0] = large-small; [1][1] = large-large;
const double Size_e[2][2] = {{2.00,2.4},{2.4,2.8}};

double typicalForce,typicalGrad,force_neta,Vcontrol;
int em=0;

double max_D,tota_D=0;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradiente conjugado
double r_k[2*N],r_k_tem[2*N];
double p_k[2*N];
double X_k[2*N];
double w[2*N];
double a_k,b_k;

double vector_vector(double A1[2*N],double B1[2*N])
{
    double h=0;
    for (int i=0; i<2*N; i++)
    {
        h=h+A1[i]*B1[i];
    }
    
    return h;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Set size of particles
void def_porte(){
    
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
// Set initial positions
void def_in_pos(){
    for (int i =0 ; i<N; i++)
    {
        rx[i]=UNI*L;
        ry[i]=UNI*L;
    }
    
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// find neighbours list
void updateNebzLists(){
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
// find neighbours list using brute force -- use just for testing
void updateNebzLists_bad(){
    
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

}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// calculate forces between particles
double fuerza_contactos; int potencial=0;
void calculateForces(){
    int i,j,m,k;
    double goo,dx,dy,r2,r,sigma;
    double temp,sigma2;
    typicalForce=0;
    
    int N_con_real=0;
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
                sigma2=(double)SizeSqrd[size[i]][size[j]];
                if ( r2 < sigma2){
                    r=sqrt(r2);
                    sigma=sqrt(sigma2);
                    
                    if (potencial==0)
                    {
                        //goo=C_e*(1-r/sigma)*(1/sigma);   // armonico
                        goo=C_e*pow((1-r/sigma),1.5)*(1/sigma); // hertzian
                    }
                    else
                    {
                        goo=C_e*(pow(sigma/r,3)-1);
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
double desplaza1(){
    
    double tem;
    max_D=0;
    double dt=1;
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
void evol_tem1(){
    for (int i =0; i<1000; i++)
    {
        max_D=desplaza1();
        tota_D=tota_D+max_D;
        if (tota_D>DR)
        {
            updateNebzLists();
            tota_D=0;
        }
        calculateForces();
        
        //printf("\n %d:  contro %0.12lf force neta %0.12lf ", i+1,Vcontrol,force_neta);
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
#define    RATIO_TOL                            1e-12

double timeStep,alpha;
int numOfStepsSinceNegativePower,itr;
double fireAdvanceTime(){
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
            if ( F_INC*timeStep < DELTA_T_MAX )
                timeStep *= F_INC;
            else
                timeStep = DELTA_T_MAX;
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
int fire(){
    
    numOfStepsSinceNegativePower = 0;
    alpha = ALPHA_START;
    timeStep = TIME_STEP;
    itr = 0;
    while ( itr < MAX_ITERATIONS && force_neta > RATIO_TOL ){
        
        max_D=fireAdvanceTime();
        tota_D=tota_D+max_D;
        if (tota_D>DR)
        {
            updateNebzLists();
            tota_D=0;
        }
        
        // printf("\n %d:  contro %0.12lf force neta %0.12lf ", itr+1,Vcontrol,force_neta);
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
void impresor(){
    FILE *fp;
    char archi[150];
    
    sprintf(archi,"result/vidrio2D_hertzian _%d_%lf_%05d.dat",N,RHO,serial);
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
     unsigned int ini, fin;
     //sscanf(argv[1],"%d",&ini);   //rho
     //sscanf(argv[2],"%d",&fin);   //rho
     //sscanf(argv[3],"%lf",&RHO);   //rho
     //L =sqrt(((1.4*1.4+1)*M_PII*DOUBLE_N)/(2*RHO));
     //HALF_L=L/2;
    
    ini=0; fin=10;
    RHO=0.925;
    L =sqrt(((1.4*1.4+1)*M_PII*DOUBLE_N)/(2*RHO));
    HALF_L=L/2;
    
    for (unsigned int tt=ini; tt<fin+ini; tt++)
    {
        potencial=0;
        tota_D=0;
        serial=tt;
        srand( serial + 123 );
        double ll=L;  printf("%lf   %d    largo caja: %g    area: %lf\n",RHO,serial,ll,ll*ll);
        
        def_porte();
        def_in_pos();
        updateNebzLists();
        calculateForces();
        evol_tem1();
        fire();    
        impresor();
    }   
}


