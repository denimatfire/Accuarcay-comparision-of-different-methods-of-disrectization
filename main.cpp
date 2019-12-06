#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

void TimeChoices(){
    cout << "MENU-for time iteration" << endl;
	cout << "1: Adam Bashfort " << endl;
	cout << "2: Runga Kutta" << endl;
	cout << "3: Exit " << endl;
	cout << "Enter your choice :";

}


struct data{
        double grid_point,x_coordinate,y_coordinate,function_old,function_new,u,v,f_n,f_nold,k1,k2,k3,k4;
        }datas[201][201];
int n,m,l,b,i,j,Schoice=0;;
double dx,dy,x_coordinate,y_coordinate,xc=-0.25,yc=0,a=5,w=1500,pi=4.0*atan(1.0),phi_max=0.0,phi_min=0.0;
void SpatialChoices(){
    cout << "MENU- for spatial schemes" << endl;
	cout << "1: FOU " << endl;
	cout << "2: CDS2" << endl;
	cout << "3: UP3 " << endl;
	cout << "4: CDS4 " << endl;
	cout << "5: UP5 " << endl;
	cout << "6: Tamm$Webb " << endl;
	cout << "7: Compact Scheme" << endl;
	cout << "8: Exit " << endl;
	cout << "Enter your choice :";
	cin>>Schoice;
	cout<<"my choice ="<<Schoice;

}
double grid_generation(){
    cout<<"grid generation started"<<"\n";
    ifstream f1;
    ofstream f2;
    f1.open("input.txt");
    f1>>n>>m>>l>>b;
    dx=double(l)/(n-1);
    dy=double(l)/(m-1);
    cout<<dx<<dy;
    for(j=0;j<m;j++){
        for (i=0;i<n;i++){
            datas[i][j].x_coordinate = -0.5+(i*dx);
            datas[i][j].y_coordinate = -0.5+(j*dy);
        }
    }
    f2.open("grid.plt");
    for(j=0;j<m;j++){
        for (i=0;i<n;i++){
            f2<<datas[i][j].x_coordinate<<"\t"<<datas[i][j].y_coordinate<<"\n";
        }
        f2<<"\n";
    }
    cout<<"grid generated"<<"\n";
}
double initial_profile(){
    cout<<"initial_profile started"<<"\n";
    ofstream f3;
    for (j=0;j<m;j++){
        for(i=0;i< n;i++){
            datas[i][j].function_old=a*exp(-w*(pow(datas[i][j].x_coordinate-xc,2)+pow(datas[i][j].y_coordinate-yc,2)));
    }}
    f3.open("function_matlab.plt");
    f3<<"zone"<<"\t"<<"i=201"<<"\t"<<"j=201"<<"\n";
    for(j=0;j<m;j++){
        for (i=0;i<n;i++)
        {
           // f3<<datas[i][j].x_coordinate<<"\t"<<datas[i][j].y_coordinate<<"\t"<<datas[i][j].function_old<<"\n";
           f3<<datas[i][j].function_old<<"\n";
        }
        f3<<"\n";
    }
    cout<<"initial_profile generated"<<"\n";
}
double velocityfield(){
    ofstream f5;
        f5.open("vfield.plt");
        //velocity
        for (j=0;j<m;j++){
            for(i=0;i<n;i++){
            datas[i][j].u = -datas[i][j].y_coordinate;
            datas[i][j].v = datas[i][j].x_coordinate;
            }
        }
    f5<<"zone"<<"\t"<<"i=201"<<"\t"<<"j=201"<<"\n";
    for(j=0;j<m;j++){
        for ( i=0;i<n;i++)
        {
            f5<<datas[i][j].x_coordinate<<"\t"<<datas[i][j].y_coordinate<<"\t"<<datas[i][j].u<<"\t"<<datas[i][j].v<<"\n";
        }
        f5<<"\n";
    }

}
double max(double x, double y){
      return x>y?x:y;
    }

double FOU(){
            cout<<"FOU called"<<"\n";
            for (j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
            // datas[i][j].f_nold=datas[i][j].f_n;
                datas[i][j].f_n = -(-max(-datas[i][j].u,0)*(datas[i+1][j].function_old - datas[i][j].function_old)/dx + max(datas[i][j].u,0)*(datas[i][j].function_old - datas[i-1][j].function_old)/dx) -(-max(-datas[i][j].v,0)*(datas[i][j+1].function_old - datas[i][j].function_old)/dy + max(datas[i][j].v,0)*(datas[i][j].function_old - datas[i][j-1].function_old)/dy);;
                }
            }
}
double CDS2(){
            cout<<"CDS2 called"<<"\n";
            for (j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
            // datas[i][j].f_nold=datas[i][j].f_n;
                datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx) + datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
                }
            }
}
double UP3(){
        cout<<"UP3 called"<<"\n";
       for(j=1;j<=(m-2);j++){
            i=1;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            i=n-2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
        }
        for(i=1;i<=(n-2);i++){
            j=1;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            j=m-2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
        }

        for (j=2;j<m-2;j++){
            for(i=2;i<n-2;i++){
           // datas[i][j].f_nold=datas[i][j].f_n;
            datas[i][j].f_n = -(-max(-datas[i][j].u,0)*((-2)*datas[i+2][j].function_old + 10*datas[i+1][j].function_old -9*datas[i][j].function_old + 2*datas[i-1][j].function_old - datas[i-2][j].function_old)/(6*dx) + max(datas[i][j].u,0)*(datas[i+2][j].function_old - 2*datas[i+1][j].function_old + 9*datas[i][j].function_old - 10*datas[i-1][j].function_old + 2*datas[i-2][j].function_old)/(6*dx)) -(-max(-datas[i][j].v,0)*((-2)*datas[i][j+2].function_old + 10*datas[i][j+1].function_old -9*datas[i][j].function_old + 2*datas[i][j-1].function_old - datas[i][j-2].function_old)/(6*dy) + max(datas[i][j].v,0)*(datas[i][j+2].function_old - 2*datas[i][j+1].function_old + 9*datas[i][j].function_old - 10*datas[i][j-1].function_old + 2*datas[i][j-2].function_old)/(6*dy));
            }
        }
}
double CDS4(){
        cout<<"CDS4 called"<<"\n";
        for(j=1;j<=(m-2);j++){
            i=1;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx) + datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            i=n-2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx) + datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
        }
        for(i=1;i<=(n-2);i++){
            j=1;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx) + datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            j=m-2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx) + datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
        }
        for(j=2;j<(m-2);j++){
            for(i=2;i<(n-2);i++){
           // datas[i][j].f_nold=datas[i][j].f_n;
            datas[i][j].f_n = -(datas[i][j].u*(-datas[i+2][j].function_old +8*datas[i+1][j].function_old - 8*datas[i-1][j].function_old + datas[i-2][j].function_old)/(12*dx) + datas[i][j].v*(-datas[i][j+2].function_old +8*datas[i][j+1].function_old - 8*datas[i][j-1].function_old + datas[i][j-2].function_old)/(12*dy));
            }
        }
}
double UP5(){
        cout<<"UP5 called"<<"\n";
       for(j=1;j<=(m-2);j++){
            i=1;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            i=2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            i=n-3;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            i=n-2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
       }
        for(i=1;i<=(n-2);i++){
            j=1;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            j=2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            j=m-3;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            j=m-2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
        }

        for (j=3;j<m-3;j++){
            for(i=3;i<n-3;i++){
           // datas[i][j].f_nold=datas[i][j].f_n;
            datas[i][j].f_n = -(-max(-datas[i][j].u,0)*( 3*datas[i+3][j].function_old -24*datas[i+2][j].function_old + 105*datas[i+1][j].function_old -20*datas[i][j].function_old - 75*datas[i-1][j].function_old + 12*datas[i-2][j].function_old - datas[i-3][j].function_old)/(120*dx) + max(datas[i][j].u,0)*(datas[i+3][j].function_old- 12*datas[i+2][j].function_old + 75*datas[i+1][j].function_old + 20*datas[i][j].function_old  - 105*datas[i-1][j].function_old + 24*datas[i-2][j].function_old - 3*datas[i-3][j].function_old)/(120*dx)) -(-max(-datas[i][j].v,0)*( 3*datas[i][j+3].function_old -24*datas[i][j+2].function_old + 105*datas[i][j+1].function_old -20*datas[i][j].function_old - 75*datas[i][j-1].function_old + 12*datas[i][j-2].function_old - datas[i][j-3].function_old)/(120*dy) + max(datas[i][j].v,0)*(datas[i][j+3].function_old- 12*datas[i][j+2].function_old + 75*datas[i][j+1].function_old + 20*datas[i][j].function_old  - 105*datas[i][j-1].function_old + 24*datas[i][j-2].function_old - 3*datas[i][j-3].function_old)/(120*dy));
            }
        }
}
double TAM(){
        cout<<"TAM called"<<"\n";
       for(j=1;j<=(m-2);j++){
            i=1;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            i=2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx))  -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            i=n-3;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            i=n-2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
        }
        for(i=1;i<=(n-2);i++){
            j=1;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            j=2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            j=m-3;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
            j=m-2;
            datas[i][j].f_n = -(datas[i][j].u*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(2*dx)) -(datas[i][j].v*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(2*dy));
        }
        for (j=3;j<m-3;j++){
            for(i=3;i<n-3;i++){
           // datas[i][j].f_nold=datas[i][j].f_n;
            datas[i][j].f_n = (-datas[i][j].u*( 0.02651995*datas[i+3][j].function_old -0.18941314*datas[i+2][j].function_old + 0.79926643*datas[i+1][j].function_old -0.79926643*datas[i-1][j].function_old + 0.18941314*datas[i-2][j].function_old - 0.02651995*datas[i-3][j].function_old) -datas[i][j].v*( 0.02651995*datas[i][j+3].function_old -0.18941314*datas[i][j+2].function_old + 0.79926643*datas[i][j+1].function_old -0.79926643*datas[i][j-1].function_old + 0.18941314*datas[i][j-2].function_old - 0.02651995*datas[i][j-3].function_old))/dx;
            }
        }
}
double COMPACT(){
        cout<<"compact called"<<"\n";
        double a[201], b[201], c[201], d[201], p[201], q[201];
        double f_primex[201][201];
        for(j=1;j<m-1;j++){
            d[0]=(datas[1][j].function_old - datas[0][j].function_old)/dx;
            d[1]=(datas[2][j].function_old - datas[0][j].function_old)/(2*dx);
            d[n-2]=(datas[n-1][j].function_old - datas[n-3][j].function_old)/(2*dx);
            d[n-1]=(datas[n-1][j].function_old - datas[n-2][j].function_old)/dx;
            a[0]=1.0;
            a[1]=1.0;
            a[n-2]=1.0;
            a[n-1]=1.0;
            b[0]=0.0;
            b[1]=0.0;
            b[n-2]=0.0;
            b[n-1]=0.0;
            c[0]=0.0;
            c[1]=0.0;
            c[n-2]=0.0;
            c[n-1]=0.0;
        for(i=2; i<n-2;i++){
            a[i]=1.0;
            b[i]= -1.0/3.0;
            c[i] = -1.0/3.0;
            d[i]=  (datas[i+2][j].function_old - datas[i-2][j].function_old)/(36*dx) + 14*(datas[i+1][j].function_old - datas[i-1][j].function_old)/(18*dx);
        }
            p[0]=b[0]/a[0];
            q[0]=d[0]/a[0];
        for(i=1; i<=n-1;i++){
            p[i]=b[i]/(a[i]-c[i]*p[i-1]);
            q[i]=(d[i]+c[i]*q[i-1])/(a[i]-c[i]*p[i-1]);
        }
        f_primex[n-1][j] = d[n-1];
        for(i=n-2;i>0;i--){
            f_primex[i][j] = q[i] + p[i]*f_primex[i+1][j];
        }
    }

double  f_primey[201][201] ;

     for(i=1;i<n-1;i++){
        d[0]=(datas[i][1].function_old - datas[i][0].function_old)/dy;
        d[1]=(datas[i][2].function_old - datas[i][0].function_old)/(2*dy);
        d[m-2]=(datas[i][m-1].function_old - datas[i][m-3].function_old)/(2*dy);
        d[m-1]=(datas[i][m-1].function_old - datas[i][m-2].function_old)/dy;
        a[0]=1.0;
        a[1]=1.0;
        a[m-2]=1.0;
        a[m-1]=1.0;
        b[0]=0.0;
        b[1]=0.0;
        b[m-2]=0.0;
        b[m-1]=0.0;
        c[0]=0.0;
        c[1]=0.0;
        c[m-2]=0.0;
        c[m-1]=0.0;
        for(j=2; j<m-2;j++){
            a[j]= 1.0;
            b[j]= -1.0/3.0;
            c[j] = -1.0/3.0;
            d[j] =  (datas[i][j+2].function_old - datas[i][j-2].function_old)/(36*dy) + 14*(datas[i][j+1].function_old - datas[i][j-1].function_old)/(18*dy);
        }
        p[0]=b[0]/a[0];
        q[0]=d[0]/a[0];
        for(j=1; j<=m-1;j++){
            p[j] = b[j]/(a[j]-c[j]*p[j-1]);
            q[j]=(d[j]+ c[j]*q[j-1])/(a[j]- c[j]*p[j-1]);
        }
        f_primey[i][m-1] = d[m-1];
        for(j=m-2;j>0;j--){
        f_primey[i][j] = q[j] + p[j]*f_primey[i][j+1];
        }
    }


        for(j=1;j<m-1;j++){
            for(i=1; i<n-1; i++){
            datas[i][j].f_n= -(datas[i][j].u*f_primex[i][j]) -(datas[i][j].v*f_primey[i][j]);
            }
        }
}
double functions(){
    cout<<"\n"<<Schoice<<"\n";
   /* if(Schoice=1){
    FOU();
    }
    else if(Schoice=2){
     CDS2();
    }
    else if(Schoice=3){
     UP3();
    }
    else if(Schoice=4){
     CDS4();
    }
    else if(Schoice=5){
     UP5();
    }
    else if(Schoice=6){
     TAM();
    }
    else if(Schoice=7){
     COMPACT();
    }*/
        switch (Schoice)
		{
		case 1:
			FOU();
			break;
		case 2:
			CDS2();
			break;
		case 3:
			UP3();
			break;
		case 4:
			CDS4();
			break;
        case 5:
			UP5();
			break;
        case 6:
			TAM();
			break;
        case 7:
			COMPACT();
			break;
		case 8:
			break;
		default:
			cout << "Invalid input" << endl;
		}
}


double AB2(){
        int count=0;
        double t=0,dt=0.0001;

        do{
            cout<<"calling function"<<"\n";
            functions();
            if (count<=100){
            cout<<"in euler explicit"<<"\n";
               for(j=1;j<m-1;j++){
                    for(i=1;i<n-1;i++){
                    datas[i][j].function_new = datas[i][j].function_old + dt*datas[i][j].f_n;
                    datas[i][j].function_old = datas[i][j].function_new;
                    }
                }
            }

            else {
                cout<<"in AB2"<<"\n";
                for(j=1;j<m-1;j++){
                    for(i=1;i<n-1;i++){
                    datas[i][j].function_new = datas[i][j].function_old + dt*(3*datas[i][j].f_n - datas[i][j].f_nold )/2;
                    datas[i][j].function_old = datas[i][j].function_new;
                    }
                 }
             }
            cout<<"no of iterartion="<<count<<" time ="<<t<<"\n";
            count++;
            t=t+dt;
            for (j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
                datas[i][j].f_nold = datas[i][j].f_n;
                }
            }

        }while(t<6*pi);

    }
double RK4(){
        int count=0;
        double t=0,dt=0.001;
        do{
            for(j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
                    datas[i][j].function_new = datas[i][j].function_old;
                }
            }
            //cout<<"calling function"<<"\n";
             functions();
            for(j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
                datas[i][j].k1 = dt*(datas[i][j].f_n );

                }
            }
                 cout<<"k1 calculated"<<"\n";
            for(j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
                datas[i][j].function_old = datas[i][j].function_new + datas[i][j].k1/2.0;
                }
            }
                functions();
            for(j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
                datas[i][j].k2 = dt*(datas[i][j].f_n);
                }
            }
                cout<<"k2 calculated"<<"\n";
            for(j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
                datas[i][j].function_old = datas[i][j].function_new + datas[i][j].k2/2.0;
                }
            }
        functions();
            for(j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
                datas[i][j].k3 = dt*(datas[i][j].f_n);
                }
            }
                cout<<"k3 calculated"<<"\n";
            for(j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
                datas[i][j].function_old = datas[i][j].function_new + datas[i][j].k3;
                }
            }
        functions();
            for(j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
                datas[i][j].k4 = dt*(datas[i][j].f_n);
                }
            }
            cout<<"k4 calculated"<<"\n";
            for(j=1;j<m-1;j++){
                for(i=1;i<n-1;i++){
                datas[i][j].function_old = datas[i][j].function_new + (datas[i][j].k1 + 2*datas[i][j].k2 + 2*datas[i][j].k3 + datas[i][j].k4)/6.0;
                //datas[i][j].function_old = datas[i][j].function_new;
                }
            }

            cout<<"no of iterartion="<<count<<" time ="<<t<<"\n";
            count++;
            t=t+dt;

        }while(t<6*pi);
}
double phi_values(){
    for(j=1;j<(m-1);j++){
        for ( i=1;i<(n-1);i++){
            if(datas[i][j].function_old>phi_max){
                    phi_max=datas[i][j].function_old;
            }
            if(datas[i][j].function_old<phi_min){
                phi_min=datas[i][j].function_old;
            }
        }
    }

    cout<<"phi max="<<phi_max<<"\t"<<"phi min="<<phi_min;
}
int main(){
    ofstream f4;
    f4.open("phi.txt");
    cout<<"grid generation starting"<<"\n";
    grid_generation();
    cout<<"initial profile starting"<<"\n";
    initial_profile();
    cout<<"going for time iteration"<<"\n";
    velocityfield();
    int Tchoice;
	do
	{
		TimeChoices();
		cin >> Tchoice;
		SpatialChoices();
		switch (Tchoice)
		{
		case 1:
			AB2();
			break;
		case 2:
			RK4();
			break;
		case 3:
			break;
		default:
			cout << "Invalid input" << endl;
		}
		 phi_values();
        for(j=0;j<m;j++){
            for ( i=0; i<n; i++){
            f4<<datas[i][j].function_old<<"\n";
            }
            f4<<"\n";
        }
	}while (Tchoice != 3);

    //f4<<"zone"<<"\t"<<"i=201"<<"\t"<<"j=201"<<"\n";

    cout<<"end";
}
