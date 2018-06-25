    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
// Author: a
#include <QcoreApplication>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

const int iter=1000000;

double force(double x, double y);

double abv(double x);

class rkn{
public:
    rkn(double x0, double y0, double px0, double py0);
    std::vector<double> psosx; //contains coordinates of trajectory in PSOS
    std::vector<double> psosp;
    void integrator();

private:
    double x[2];
    double p[2];
};

rkn::rkn(double x0, double y0, double px0, double py0)
{
    x[0]=x0;
    x[1]=y0;
    p[0]=px0;
    p[1]=py0;
}

void rkn::integrator()
{
    double c[4]={0,1/5.0,2/3.0,1};
    double a[4][4]={{0,0,0,0},{1/50.0,0,0,0},{-1/27.0,7/27.0,0,0},{3/10.0,-2/35.0,9/35.0,0}};
    double b[4]={14/336.0,100/336.0,54/336.0,0};
    double d[4]={14/336.0,25/336.0,162/336.0,35/336.0}; //RKN coefficients
    int current=0; //counter for the PSOS data
    double k[4][2];
    double sk[2]; //just for convenience
    double dt=0.005; //integration stepsize
    for(int n=0;n<=iter;n++)
    {
        if(abv(x[1])<=0.001) //check if trajectory crosses PSOS
        {
            //std::cout<<abv(x[1])<<"\n";
            psosx.insert(psosx.end(),x[0]);
            psosp.insert(psosp.end(),p[0]);
            current+=1;
        }
        for(int i=0;i<=3;i++)
        {
            k[i][0]=0;
            k[i][1]=0;
            sk[0]=0;
            sk[1]=0;
            if(i!=0)
            {
                for(int j=0;j<=i-1;j++)
                {
                    sk[0]+=k[j][0]*a[i][j];
                    sk[1]+=k[j][1]*a[i][j];
                }
            }
            k[i][0]=force(x[0]+p[0]*dt*c[i]+sk[0]*dt*dt,x[1]+p[1]*dt*c[i]+sk[1]*dt*dt);
            k[i][1]=force(x[1]+p[1]*dt*c[i]+sk[0]*dt*dt,x[0]+p[0]*dt*c[i]+sk[1]*dt*dt); //force for rkn
        }
        for(int i=0;i<=1;i++)
        {
            x[i]+=p[i]*dt+(k[0][i]*b[0]+k[1][i]*b[1]+k[2][i]*b[2])*dt*dt;
            p[i]+=(k[0][i]*d[0]+k[1][i]*d[1]+k[2][i]*d[2]+k[3][i]*d[3]+k[4][i]*d[4])*dt; //integration via rkn
        }
    }
}

double force(double x,double y)
{
    return -x+y*y*y-3*x*y*y; //from the equation of motion
}

double abv(double x)
{
    if(x>=0)
        return x;
    else
        return -x;
}

/*Vector2f force(Vector2f x) //some strange issues make using the Vector2f class not viable
{
    Vector2f f;
    f(0)=-x(0)-3*x(0)*x(1)*x(1)+x(0)*x(0)*x(0);
    f(1)=-x(1)-3*x(0)*x(0)*x(1)+x(1)*x(1)*x(1);
    return f;
}*/

int main()
{
    rkn system(2,0,0,1);
    system.integrator();
    std::ofstream ofs;
    ofs.open("output.txt"); //output
    for(unsigned int i=0;i<=system.psosx.size()-1;i++)
    {
        ofs<<system.psosx[i]<<"\t"<<system.psosp[i]<<"\n";
    }
    std::cout<<"done";
}
