#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

//Global Variables
int i,j,I,J,nh=126,nv=126,n=nv+2;
double dx=1.0/n;
double dy=1.0/n;
double dt = 0.0001;
double Pr = 0.71;
double Ra = 10000;
double S = 0.0;
double error_p_prime =0.0;
double limit = 1e-8;
double t = 0.0;
int temp_counter = 0.0;
int count = 0.0;
double error =0.0;
float w=1.0;
double a_prime;
float w1 =1.0;
int count_time = 0;
int count_inner_loop =0;
float Re = 100.0;

double p=dt/(2*dx*dx),q=dt/(2*dy*dy);


//Structure
struct data{
        double P_x,P_y,U_x,U_y,V_x,V_y;
        double U,V,P,T,U_old,V_old;
        double ue2,uw2,du2dx,un,vn,us,vs,duvdy,Cu;
        double ue2_old,uw2_old,du2dx_old,un_old,vn_old,us_old,vs_old,duvdy_old,Cu_old;
        double d2udx2,d2udy2,Du;
        double d2udx2_old,d2udy2_old,Du_old;
        double vn2,vs2,dv2dy,ue,ve,uw,vw,duvdx,Cv;
        double vn2_old,vs2_old,dv2dy_old,ue_old,ve_old,uw_old,vw_old,duvdx_old,Cv_old;
        double d2vdx2,d2vdy2,Dv;
        double d2vdx2_old,d2vdy2_old,Dv_old;
        double U_tilda,V_tilda,P_tilda;//,T_prev,T_new;
        //double Ap,As,Aw,Ae,An,Bp;
        //double Ap_prev,As_prev,Aw_prev,Ae_prev,An_prev,Bp_prev;
        //double ap,ae,aw,an,as;
        double rho_temp,rho;
        double AP,AS,AW,AN,AE,BP,P_tilda_new,RHS;
        }node[200][200];

//declearing all subroutines
void printer();
//void grig_generation();
//void boundary_condition_energy_equation2();
//void ab2_cn_energy_equation_coefficients();
//void explicit_energy_equation_solver();
//void boundary_condition_energy_equation1();
//void coefficients_energy_equation1();
void equation_3b();
void equation_3a();
void calculate_RHS_equation_6();
void Velocity_BC();
void equation_1a();
void equation_1b();
void Y_momentum_solver();
void X_momentum_solver();
void intialisation();
void grid_output();
void MAC();
void correct_pressure();
void pressure_solver();
void pressure_coffecients();
void swap();
void store_old_velocities();
void X_momentum_solver_old();
void Y_momentum_solver_old();
void equation_1a_ab2();
void equation_1b_ab2();

//main function
int main()
{
   // grig_generation();
    intialisation();
   // grid_output();
    MAC();
}

//

//Print grid in file


//intialising variables
void intialisation(){

    for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
            node[i][j].U = 0.0;//U initialised
			node[i][j].U_old = 0.0;//U initialised
            node[i][j].V = 0.0;// V initialised
			node[i][j].V_old = 0.0;// V initialised
            node[i][j].P = 1.0;// P initialised
           // node[i][j].P_tilda = 0.0;//pressure correction
            node[i][j].T = 0.0;//Temp inilialised
           // node[i][j].T_prev = 0.0;//n-1 time level temperature
        }
    }

}

void MAC(){

        

		do{
            count_inner_loop =0;
        	count_time++;
        	
          		
          		
          	
          		Velocity_BC();
          		store_old_velocities();
			equation_1a();//pedicting u velocity
			equation_1b();//predicting v velocity
		
			//Velocity_BC();
		 	//equation_1a_ab2();
		 	//equation_1b_ab2();
		 	//store_old_velocities();
			calculate_RHS_equation_6();//calculating pressure correction
			pressure_coffecients();
			

				do{
					pressure_solver();
					equation_3a();//correcting u velocity
					equation_3b();//correcting v velocity
					calculate_RHS_equation_6();
                        error_p_prime = sqrt(S/(n*n));
                        count_inner_loop++;
                        cout<<"error = "<<error_p_prime<<"\n"<<" time "<<count_inner_loop<<"\n";

						//break ;
				}while(error_p_prime>limit);
				swap();
           // correct_pressure();

			
			t=t+dt;

			cout<<"time = "<<count_time<<"\n";
		}while(t<=1);
					printer();

}

void X_momentum_solver(){
    //X momentum solved in ucell
    
    // convective term
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].ue2 =pow(node[i+1][j].U + node[i][j].U,2);
                node[i][j].uw2 = pow(node[i-1][j].U + node[i][j].U,2);
                node[i][j].du2dx = (node[i][j].ue2 - node[i][j].uw2)/(4*dx);

                node[i][j].un = (node[i][j+1].U + node[i][j].U);
                node[i][j].us = (node[i][j-1].U + node[i][j].U);
                node[i][j].vn = (node[i+1][j].V + node[i][j].V);
                node[i][j].vs = (node[i][j-1].V + node[i+1][j-1].V);
                node[i][j].duvdy = (node[i][j].un*node[i][j].vn - node[i][j].us*node[i][j].vs)/(4*dy);

                node[i][j].Cu = node[i][j].du2dx + node[i][j].duvdy;

        }
    }
    //Diffusive term
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].d2udx2 = (node[i+1][j].U - 2.0*node[i][j].U + node[i-1][j].U)/(dx*dx);
                node[i][j].d2udy2 = (node[i][j+1].U - 2.0*node[i][j].U + node[i][j-1].U)/(dy*dy);

                node[i][j].Du = node[i][j].d2udx2 + node[i][j].d2udy2;
        }
    }
   
}

void Y_momentum_solver(){
    //Y momentum solved in V cell
  
    // convective term
    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
            node[i][j].vn2 = pow(node[i][j+1].V + node[i][j].V,2);
            node[i][j].vs2 = pow(node[i][j-1].V + node[i][j].V,2);
            node[i][j].dv2dy = (node[i][j].vn2 -node[i][j].vs2)/(4*dy);

            node[i][j].ue = (node[i][j].U + node[i][j+1].U);
            node[i][j].uw = (node[i-1][j].U + node[i-1][j+1].U);
            node[i][j].ve = (node[i][j].V + node[i+1][j].V);
            node[i][j].vw = (node[i-1][j].V + node[i][j].V);
            node[i][j].duvdx = (node[i][j].ue*node[i][j].ve - node[i][j].uw*node[i][j].vw)/(4*dx);

            node[i][j].Cv = node[i][j].dv2dy + node[i][j].duvdx;
        }
    }
    //diffusive term
    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
                node[i][j].d2vdx2 = (node[i+1][j].V - 2.0*node[i][j].V + node[i-1][j].V)/(dx*dx);
                node[i][j].d2vdy2 = (node[i][j+1].V - 2.0*node[i][j].V + node[i][j-1].V)/(dy*dy);

                node[i][j].Dv = node[i][j].d2vdx2 + node[i][j].d2vdy2;
        }
    }
  

}

//predictor step

//U momentum in U cell
void equation_1a(){

	X_momentum_solver();
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].U_tilda = node[i][j].U - dt*node[i][j].Cu + (1.0/Re)*dt*node[i][j].Du - (node[i+1][j].P - node[i][j].P)*(dt/dx);
        }
    }
  
}
//V momentum in V cell

void equation_1b(){
   
	Y_momentum_solver();
    for(j=1;j<=n-1;j++){
            for(i=1;i<=n;i++){
                    node[i][j].V_tilda = node[i][j].V - dt*node[i][j].Cv + (1.0/Re)*dt*node[i][j].Dv - (node[i][j+1].P - node[i][j].P)*(dt/dy);// + dt*Ra*Pr*(node[i][j].T);
            }
    }
   
}
//boundary conditions of velocity
void Velocity_BC(){

  
    //west wall
    for(j=0;j<=n;j++){
        node[0][j].U = 0.0;//no slip
		node[0][j].U_old = 0.0;//no slip
        node[0][j].U_tilda = 0.0;
        node[0][j].V = -node[1][j].V;
		node[0][j].V_old = -node[1][j].V_old;
        node[0][j].V_tilda = -node[1][j].V_tilda;//no slip
    }
    //east wall
     for(j=0;j<=n;j++){
        node[n][j].U = 0.0;//no slip
		node[n][j].U_old = 0.0;//no slip
        node[n][j].U_tilda = 0.0;//no slip
        node[n+1][j].V = -node[n][j].V;//no slip
		node[n+1][j].V_old = -node[n][j].V_old;//no slip
        node[n+1][j].V_tilda = -node[n][j].V_tilda;
     }
     //top wall
     for(i=0;i<=n;i++){
        node[i][n+1].U = 2*1.0 - node[i][n].U;//no slip
		node[i][n+1].U_old =2*1.0 - node[i][n].U_old;//no slip
        node[i][n+1].U_tilda =2*1.0 - node[i][n].U_tilda;
        node[i][n].V = 0.0;//no slip
		node[i][n].V_old = 0.0;//no slip
        node[i][n].V_tilda = 0.0;//no slip
     }
     //bottom wall
     for(i=0;i<=n;i++){
        node[i][0].U = - node[i][1].U;//no slip
		node[i][0].U_old = - node[i][1].U;//no slip
        node[i][0].U_tilda = - node[i][1].U_tilda;//no slip
        node[i][0].V = 0.0;//no slip
		node[i][0].V_old = 0.0;//no slip
        node[i][0].V_tilda = 0.0;//no slip
     }
   
}
//calculation of pressure correction
void calculate_RHS_equation_6(){
    
                a_prime = -2.0*dt*(1.0/pow(dx,2) + 1.0/pow(dy,2) );
              
                S = 0.0;
    for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
                node[i][j].RHS = w1*((node[i][j].U_tilda - node[i-1][j].U_tilda)/dx + (node[i][j].V_tilda - node[i][j-1].V_tilda)/dy)/a_prime;
                node[i][j].BP = (1.0/dt)*(node[i][j].U_tilda - node[i-1][j].U_tilda)/dx + (node[i][j].V_tilda - node[i][j-1].V_tilda)/dy;
                S = S + node[i][j].RHS;
        }
    }
	
	
  
}
// u correction in u cell
void equation_3a(){
   
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].U_tilda = node[i][j].U_tilda - (dt/dx)*(node[i+1][j].P_tilda - node[i][j].P_tilda);
        }
    }
   
}
// v correction in v cell
void equation_3b(){
   
    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
                node[i][j].V_tilda = node[i][j].V_tilda - (dt/dy)*(node[i][j+1].P_tilda - node[i][j].P_tilda);
        }
    }
   
}


//subroutin to check modules
void printer(){
ofstream f2;
f2.open("check.txt");
for(j=0;j<=n+1;j++){
        for(i=0;i<=n+1;i++){
     f2<<i<<"\t"<<j<<"\t"<<node[i][j].U<<"\t"<<node[i][j].V<<"\t"<<node[i][j].P_tilda<<"\n";
        }
}
}
void correct_pressure(){
   
	for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
			node[i][j].P = node[i][j].P + node[i][j].P_tilda;
		}
	}
}


void pressure_coffecients(){
    for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
            node[i][j].AS = 0.0;
            node[i][j].AN = 0.0;
            node[i][j].AE = 0.0;
            node[i][j].AW = 0.0;
        }
    }
    for(j=2;j<=n;j++){
        for(i=1;i<=n;i++){
            node[i][j].AS = 1.0/(dy*dy);
        }
    }

    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
            node[i][j].AN = 1.0/(dy*dy);
        }
    }

    for(j=1;j<=n;j++){
        for(i=2;i<=n;i++){
            node[i][j].AW = 1.0/(dx*dx);
        }
    }

    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
            node[i][j].AE = 1.0/(dy*dy);
        }
    }
    for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
            node[i][j].AP = -(node[i][j].AW + node[i][j].AS + node[i][j].AN + node[i][j].AE);
        }
    }
}

void pressure_solver(){
    for(j=1;j<=n;j++){
        for(i=1;i<=n;i++){
            node[i][j].rho = node[i][j].BP - (node[i][j].AW*node[i-1][j].P_tilda + node[i][j].AE*node[i+1][j].P_tilda + node[i][j].AN*node[i][j+1].P_tilda + node[i][j].AS*node[i][j-1].P_tilda + node[i][j].AP*node[i][j].P_tilda);
            node[i][j].P_tilda_new = node[i][j].P_tilda + w*node[i][j].rho/node[i][j].AP;
                //s = s+ pow((node[i][j].P_tilda_new-node[i][j].P_tilda), 2);
                node[i][j].P_tilda = node[i][j].P_tilda_new;
        }
    }
}
void swap(){
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].U = node[i][j].U_tilda;
        }
    }
	    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
                node[i][j].V = node[i][j].V_tilda ;
        }
    }
}
void store_old_velocities(){
	for(j=1;j<=n;j++){
		for(i=1;i<=n-1;i++){
		node[i][j].U_old = node[i][j].U;
	
		}
	}
	
	for(j=1;j<=n-1;j++){
		for(i=1;i<=n;i++){
		node[i][j].V_old = node[i][j].V;
		}
	}
}

void X_momentum_solver_old(){
    //X momentum solved in ucell
    
    // convective term
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].ue2_old =pow((node[i+1][j].U_old + node[i][j].U_old),2);
                node[i][j].uw2_old = pow((node[i-1][j].U_old + node[i][j].U_old),2);
                node[i][j].du2dx_old = (node[i][j].ue2_old - node[i][j].uw2_old)/(4*dx);

                node[i][j].un_old = (node[i][j+1].U_old + node[i][j].U_old);
                node[i][j].us_old = (node[i][j-1].U_old + node[i][j].U_old);
                node[i][j].vn_old = (node[i+1][j].V_old + node[i][j].V_old);
                node[i][j].vs_old = (node[i][j-1].V_old + node[i+1][j-1].V_old);
                node[i][j].duvdy_old = (node[i][j].un_old*node[i][j].vn_old - node[i][j].us_old*node[i][j].vs_old)/(4*dy);

                node[i][j].Cu_old = node[i][j].du2dx_old + node[i][j].duvdy_old;

        }
    }
    //Diffusive term
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].d2udx2_old = (node[i+1][j].U_old - 2.0*node[i][j].U_old + node[i-1][j].U_old)/(dx*dx);
                node[i][j].d2udy2_old = (node[i][j+1].U_old - 2.0*node[i][j].U_old + node[i][j-1].U_old)/(dy*dy);

                node[i][j].Du_old = node[i][j].d2udx2_old + node[i][j].d2udy2_old;
        }
    }
   
}
void Y_momentum_solver_old(){
    //Y momentum solved in V cell
  
    // convective term
    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
            node[i][j].vn2_old = pow((node[i][j+1].V_old + node[i][j].V_old),2);
            node[i][j].vs2_old = pow((node[i][j-1].V_old + node[i][j].V_old),2);
            node[i][j].dv2dy_old = (node[i][j].vn2_old -node[i][j].vs2_old)/(4*dy);

            node[i][j].ue_old = (node[i][j].U_old + node[i][j+1].U_old);
            node[i][j].uw_old = (node[i-1][j].U_old + node[i-1][j+1].U_old);
            node[i][j].ve_old = (node[i][j].V_old + node[i+1][j].V_old);
            node[i][j].vw_old = (node[i-1][j].V_old + node[i][j].V_old);
            node[i][j].duvdx_old = (node[i][j].ue_old*node[i][j].ve_old - node[i][j].uw_old*node[i][j].vw_old)/(4*dx);

            node[i][j].Cv_old = node[i][j].dv2dy_old + node[i][j].duvdx_old;
        }
    }
    //diffusive term
    for(j=1;j<=n-1;j++){
        for(i=1;i<=n;i++){
                node[i][j].d2vdx2_old = (node[i+1][j].V_old - 2.0*node[i][j].V_old + node[i-1][j].V_old)/(dx*dx);
                node[i][j].d2vdy2_old = (node[i][j+1].V_old - 2.0*node[i][j].V_old + node[i][j-1].V_old)/(dy*dy);

                node[i][j].Dv_old = node[i][j].d2vdx2_old + node[i][j].d2vdy2_old;
        }
    }
}

void equation_1a_ab2(){

	X_momentum_solver();
	X_momentum_solver_old();
    for(j=1;j<=n;j++){
        for(i=1;i<=n-1;i++){
                node[i][j].U_tilda = node[i][j].U - (dt/2.0)*(3.0*node[i][j].Cu - node[i][j].Cu_old) + (1.0/Re)*(dt/2.0)*(3.0*node[i][j].Du - node[i][j].Du_old);
        }
    }
  
}

void equation_1b_ab2(){
   
	Y_momentum_solver();
	Y_momentum_solver_old();
    for(j=1;j<=n-1;j++){
            for(i=1;i<=n;i++){
                    node[i][j].V_tilda = node[i][j].V - (dt/2.0)*(3.0*node[i][j].Cv - node[i][j].Cv_old) + (1.0/Re)*(dt/2.0)*(3.0*node[i][j].Dv - node[i][j].Dv_old); //(node[i][j+1].P - node[i][j].P)*(dt/dy);// + dt*Ra*Pr*(node[i][j].T);
            }
    }
   
}

