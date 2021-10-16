#include<stdio.h>
#include<math.h>

#include"soft1.h"

double grid();
double calculation();
double print();
double vorticity_calculation();
double vorticity_boundary_condition();
double guess();
double velocity_calculation();
double velocity_boundary_condition();
double error_calculation();
double update();
double stream_calculation();
double stream_boundary_condition();

double applied_field_initial();
double applied_field_boundary_condition();
double applied_field_calculation();
double applied_field_error();
double applied_field_update();
double applied_field_print();
double applied_field();

double induced_field_initial();
double induced_field_boundary_condition();
double induced_field_calculation();
double induced_field_error();
double induced_field_update();
double induced_field_print();
double induced_field();

double flow_field();
double source_calculation();
double print_source();

double species_initial();
double species_boundary_condition();
double species_calculation();
double species_error();
double species_update();
double species_print();
double species_transport();

//*************************************************************************************************************************************************************************//
																	//MAIN FUNCTION//
//*************************************************************************************************************************************************************************//
void main()
{
	grid();
	calculation();
	
	m1 = 1.3*m/5;
	m2 = 2.5*m/5;
	m3 = 3.7*m/5;
	m4 = 0.5*m/5;
	m5 = 4.5*m/5;
	
	
	FILE *fp1;
	fp1 = fopen("cc16.dat","w");
	fprintf(fp1,"\n%f\t\t%f",u[m1][n2],sp[m1][n2]);
	fprintf(fp1,"\n%f\t\t%f",u[m2][n2],sp[m2][n2]);
	fprintf(fp1,"\n%f\t\t%f",u[m3][n2],sp[m3][n2]);
	fprintf(fp1,"\n%f\t\t%f",u[m4][n/2],sp[m4][n/2]);
	fprintf(fp1,"\n%f\t\t%f",u[m1][n/2],sp[m1][n/2]);
	fprintf(fp1,"\n%f\t\t%f",u[m2][n/2],sp[m2][n/2]);
	fprintf(fp1,"\n%f\t\t%f",u[m3][n/2],sp[m3][n/2]);
	fprintf(fp1,"\n%f\t\t%f",u[m5][n/2],sp[m5][n/2]);
	fprintf(fp1,"\n%f\t\t%f",u[m1][n1],sp[m1][n1]);
	fprintf(fp1,"\n%f\t\t%f",u[m2][n1],sp[m2][n1]);
	fprintf(fp1,"\n%f\t\t%f",u[m3][n1],sp[m3][n1]);
	
	fclose(fp1);
}
//*************************************************************************************************************************************************************************//
																	//TOTAL CALCULATION//
//*************************************************************************************************************************************************************************//
double calculation()
{
	app_pot = length;
	printf("\ncalculating the induced field...");
	induced_field();	
	/*printf("\ncalculating the applied field...");
	applied_field();*/
	source_calculation();
	print_source();
	printf("\ncalculating the flow field...");
	flow_field();
	printf("\ncalculating the species transport...");
	species_transport();
	print();
}
//*************************************************************************************************************************************************************************//
																	//PRINTING VALUES//
//*************************************************************************************************************************************************************************//
double print()
{
	FILE *fp1;
	FILE *fp2;
	FILE *fp3;
	FILE *fp4;
	fp1 = fopen("vorticity.dat","w");
	fp2 = fopen("stream.dat","w");
	fp3 = fopen("u.dat","w");
	fp4 = fopen("v.dat","w");
	fprintf(fp1,"ZONE\tJ=%d\tI=%d",m+1,n+1);
	fprintf(fp2,"ZONE\tJ=%d\tI=%d",m+1,n+1);
	fprintf(fp3,"ZONE\tJ=%d\tI=%d",m+1,n+1);
	fprintf(fp4,"ZONE\tJ=%d\tI=%d",m+1,n+1);
	for(i=0;i<=m;i++)
	{
		for(j=0;j<=n;j++)
		{
			fprintf(fp1,"\n%f\t%f\t%f",x[i][j],y[i][j],omega[i][j]);
			fprintf(fp2,"\n%f\t%f\t%f",x[i][j],y[i][j],stream[i][j]);
			fprintf(fp3,"\n%f\t%f\t%f\t%f",x[i][j],y[i][j],u[i][j],v[i][j]);
			fprintf(fp4,"\n%f\t%f\t%f",x[i][j],y[i][j],v[i][j]);
		}
	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	
	
	//1-D Printing Code//
	
	FILE *fp5;
	FILE *fp6;
	fp5 = fopen("v500x100.dat","w");
	fp6 = fopen("c500x100.dat","w");
	int p = 4*m/5; //For velocity and species- 4th position
	for(j=0;j<=n;j++)
	{
		fprintf(fp5,"\n%f\t%f",y[0][j],u[p][j]);
		fprintf(fp6,"\n%f\t%f",y[0][j],sp[p][j]);
	}
	fclose(fp5);
	fclose(fp6);
	
}
//*************************************************************************************************************************************************************************//
																	//GRID GENERATIONS//
//*************************************************************************************************************************************************************************//
double grid()
{
	dx = (length)/m;
	dy = (height)/n;
	n1 = d*(n/2) - 1;
	n2 = n - n1;
	
	for(i=0;i<=m;i++)
	{
		for(j=0;j<=n;j++)
		{
			x[i][j] = (length)*i/m;
			y[i][j] = (height*j/n) - (height/2.0);
		}
	}
}
//*************************************************************************************************************************************************************************//
															//CALCULATING BOUNDARY CONDITIONS//
//*************************************************************************************************************************************************************************//
double vorticity_boundary_condition()
{
	//for inlet//
	for(j=0;j<=n;j++)
	{
		omega[0][j] = 2.0*(stream[0][j] - stream[1][j])/(dx*dx);
		omegaold[0][j] = 2.0*(stream[0][j] - stream[1][j])/(dx*dx);
	}
	// for outlet//
	for(j=0;j<=n;j++)
	{
		omega[m][j] = omega[m-1][j];
		omegaold[m][j] = omegaold[m-1][j];
	}
	//for upper wall//
	for(i=0;i<=m;i++)
	{
		omega[i][n] = 2.0*(stream[i][n] - stream[i][n-1])/(dy*dy);
		omegaold[i][n] = 2.0*(stream[i][n] - stream[i][n-1])/(dy*dy);
	}
	//for lower wall//
	for(i=0;i<=m;i++)
	{
		omega[i][0] = 2.0*(stream[i][0] - stream[i][1])/(dy*dy);
		omegaold[i][0] = 2.0*(stream[i][0] - stream[i][1])/(dy*dy);
	}
}

double velocity_boundary_condition()
{
	//for upper and lower wall//
	for(i=0;i<=m;i++)
	{
		u[i][n] = 0.0;
		v[i][n] = 0.0;
		
		u[i][0] = 0.0;
		v[i][0] = 0.0;
	}
	//for inlet and outlet//
	for(j=0;j<=n;j++)
	{
		u[0][j] = uin;
		v[0][j] = 0.0;
		
		u[m][j] = u[m-1][j];
		v[m][j] = v[m-1][j];
	}
}

double stream_boundary_condition()
{
	//for upper wall//
	for(i=0;i<=m;i++)
	{
		stream[i][n] = uin;
	}
	//for lower wall//
	for(i=0;i<=m;i++)
	{
		stream[i][0] = 0.0;
	}
	//for inlet//
	for(j=1;j<=n;j++)
	{
		stream[0][0] = 0.0;
		stream[0][j] = stream[0][j-1] + (uin*dy);
	}
	//for outlet//
	for(j=1;j<=n;j++)
	{
		stream[0][0] = 0.0;
		stream[m][j] = (2.0*stream[m-1][j]) - stream[m-2][j];
	}
}
//*************************************************************************************************************************************************************************//
																	//GUESSING THE VALUES//
//*************************************************************************************************************************************************************************//
double guess()
{
	for(i=0;i<=m;i++)
	{
		for(j=0;j<=n;j++)
		{
			omega[i][j] = 0.0;
			omegaold[i][j] = 0.0;
			
			u[i][j] = 0.0;
			v[i][j] = 0.0;
			stream[i][j] = 0.0;
		}
	}
}
//*************************************************************************************************************************************************************************//
															//CALCULATING VELOCITY, STREAM FUNCTION AND VORTICITY//
//*************************************************************************************************************************************************************************//
double velocity_calculation()
{
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			u[i][j] = (stream[i][j+1] - stream[i][j-1])/(2.0*dy);
			v[i][j] = (-stream[i+1][j] + stream[i-1][j])/(2.0*dx);
		}
	}
}

double stream_calculation()
{
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			stream[i][j] = ( ((stream[i+1][j] + stream[i-1][j])/(dx*dx)) + ((stream[i][j+1] + stream[i][j-1])/(dy*dy)) + omega[i][j] )/( (2.0/(dx*dx)) + (2.0/(dy*dy)) );
		}
	}
}

double vorticity_calculation()
{
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			if(j>n1&&j<n2)
			{
				flag = 0.0;
			}
			else
			{
				flag = 1.0;
			}
			fakeomega[i][j] = ( -((u[i][j]*re/(2.0*dx))*(omegaold[i+1][j] - omega[i-1][j])) - ((v[i][j]*re/(2.0*dy))*(omegaold[i][j+1] - omega[i][j-1])) 
			+ ((1.0/(dx*dx))*(omegaold[i+1][j] + omega[i-1][j])) + ((1.0/(dy*dy))*(omegaold[i][j+1] + omega[i][j-1])) 
			- ((1.0/(2.0*dy))*(f[i][j+1] - f[i][j-1])) )/( (2.0/(dx*dx)) + (2.0/(dy*dy)) + (flag*alpha*alpha) );
			omega[i][j] = omegaold[i][j] + (womega*(fakeomega[i][j] - omegaold[i][j]));
		}		
	}
}
//*************************************************************************************************************************************************************************//
																		//CALCULATING ERROR//
//*************************************************************************************************************************************************************************//
double error_calculation()
{
	sum = 0.0;
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			sum = sum + ((omega[i][j] - omegaold[i][j])*(omega[i][j] - omegaold[i][j])) ;
		}
	}
	error = sqrt(sum/((m-1)*(n-1)));
}
//*************************************************************************************************************************************************************************//
																	//UPDATING THE VALUES//
//*************************************************************************************************************************************************************************//
double update()
{
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			omegaold[i][j] = omega[i][j];
		}
	}
}
//*************************************************************************************************************************************************************************//
																		//CALCULATING BODY FORCE//
//*************************************************************************************************************************************************************************//
double source_calculation()
{
	
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			f[i][j] = -kappa*kappa*sinh(sy[i][j])*(app_pot/length);
		}
	}
	for(j=1;j<n;j++)
	{
		f[0][j] = -kappa*kappa*sinh(sy[0][j])*(app_pot/length);
		f[m][j] = -kappa*kappa*sinh(sy[m][j])*(app_pot/length);
	}
	for(i=1;i<m;i++)
	{
		f[i][0] = -kappa*kappa*sinh(sy[i][0])*(app_pot/length);
		f[i][n] = -kappa*kappa*sinh(sy[i][n])*(app_pot/length);
	}
}

double print_source()
{
	FILE *fp6;
	fp6 = fopen("source.dat","w");
	fprintf(fp6,"zone\tj=%d\ti=%d",m+1,n+1);
	for(i=0;i<=m;i++)
	{
		for(j=0;j<=n;j++)
		{
			fprintf(fp6,"\n%f\t\t%f\t\t%f",x[i][j],y[i][j],f[i][j]);
		}
	}
	fclose(fp6);
}
//*************************************************************************************************************************************************************************//
																		//CALCULATING APPLIED FIELD//
//*************************************************************************************************************************************************************************//
double applied_field()
{
	applied_field_initial();
	app_error = 1000.0;
	count = 1;
	while(app_error>1.0e-10)
	{
		applied_field_boundary_condition();
		applied_field_calculation();
		applied_field_error();
		if(count%1000==0.0)
		{
			printf("\n%d\t\t%.15f",count,app_error);
		}
		applied_field_update();
		count++;
	}
	applied_field_print();
}

double applied_field_initial()
{
	for(i=0;i<=m;i++)
	{
		for(j=0;j<=n;j++)
		{
			phi[i][j] = 0.0;
			phiold[i][j] = 0.0;
		}
	}
}

double applied_field_boundary_condition()
{
	//at lower and upper wall//
	for(i=0;i<=m;i++)
	{
		phi[i][0] = phi[i][1];
		phiold[i][0] = phiold[i][1];
		phi[i][n] = phi[i][n-1];
		phiold[i][n] = phiold[i][n-1];
	}
	//at left and right wall//
	for(j=0;j<=n;j++)
	{
		phi[0][j] = 0.0;
		phiold[0][j] = 0.0;
		phi[m][j] = app_pot;
		phiold[m][j] = app_pot;
	}
}

double applied_field_calculation()
{
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			fakephi[i][j] = (((phiold[i+1][j] + phi[i-1][j])/(dx*dx)) + ((phiold[i][j+1] + phi[i][j-1])/(dy*dy)))/(2.0*((1.0/(dx*dx)) + (1.0/(dy*dy))));
			phi[i][j] = phiold[i][j] + (wphi*(fakephi[i][j] - phiold[i][j]));
		}
	}
}

double applied_field_error()
{
	sum = 0.0;
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			sum = sum + ((phi[i][j] - phiold[i][j])*(phi[i][j] - phiold[i][j]));
		}
	}
	app_error = sqrt(sum/(m*n));
}

double applied_field_update()
{
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			phiold[i][j] = phi[i][j];
		}
	}
}

double applied_field_print()
{
	FILE *fp6;
	fp6 = fopen("applied field.dat","w");
	fprintf(fp6,"zone\tj=%d\ti=%d",m+1,n+1);
	for(i=0;i<=m;i++)
	{
		for(j=0;j<=n;j++)
		{
			fprintf(fp6,"\n%f\t\t%f\t\t%f",x[i][j],y[i][j],phi[i][j]);
		}
	}
	fclose(fp6);
}

//*************************************************************************************************************************************************************************//
																		//CALCULATING INDUCED FIELD//
//*************************************************************************************************************************************************************************//
double induced_field()
{
	app_error = 1000.0;
	count = 1;
	induced_field_initial();
	while(app_error>1.0e-15)
	{
		induced_field_boundary_condition();
		induced_field_calculation();
		induced_field_error();
		if(count%1000==0.0)
		{
			printf("\n%d\t\t%.20f",count,app_error);
		}
		count++;
		induced_field_update();
	}
	induced_field_print();
}

double induced_field_initial()
{
	for(i=0;i<=m;i++)
	{
		for(j=0;j<=n;j++)
		{
			sy[i][j] = 0.0;
			syold[i][j] = 0.0;
		}
	}
}

double induced_field_boundary_condition()
{
	
	m1 = m/5;
	m2 = 8*m/25;
	m3 = 11*m/25;
	m4 = 14*m/25;
	m5 = 17*m/25;
	m6 = 20*m/25;
	
	//at upper wall and lower wall//
	for(i=0;i<m1;i++)
	{
		sy[i][n] = 0.0;
		sy[i][0] = 0.0;
		syold[i][n] =0.0;
		syold[i][0] = 0.0;
	}
		
	for(i=m1;i<=m2;i++)
	{
		sy[i][n] = ((4.0*sy[i][n-1]) - sy[i][n-2])/3.0;
		sy[i][0] = 0.0;
		syold[i][n] = ((4.0*syold[i][n-1]) - syold[i][n-2])/3.0;
		syold[i][0] = 0.0;	
	}
	
	for(i=m2+1;i<m3;i++)
	{
		sy[i][n] = 0.0;
		sy[i][0] = ((4.0*sy[i][1]) - sy[i][2])/3.0;
		syold[i][n] =0.0;
		syold[i][0] = ((4.0*syold[i][1]) - syold[i][2])/3.0;
	}
		
	for(i=m3;i<=m4;i++)
	{
		sy[i][n] = ((4.0*sy[i][n-1]) - sy[i][n-2])/3.0;
		sy[i][0] = 0.0;
		syold[i][n] = ((4.0*syold[i][n-1]) - syold[i][n-2])/3.0;
		syold[i][0] = 0.0;	
	}
	
	for(i=m4+1;i<m5;i++)
	{
		sy[i][n] = 0.0;
		sy[i][0] =((4.0*sy[i][1]) - sy[i][2])/3.0;
		syold[i][n] =0.0;
		syold[i][0] = ((4.0*syold[i][1]) - syold[i][2])/3.0;
	}
	
	for(i=m5;i<=m6;i++)
	{
		sy[i][n] = ((4.0*sy[i][n-1]) - sy[i][n-2])/3.0;
		sy[i][0] = 0.0;
		syold[i][n] = ((4.0*syold[i][n-1]) - syold[i][n-2])/3.0;
		syold[i][0] = 0.0;	
	}
	
	for(i=m6+1;i<=m;i++)
	{
		sy[i][n] = 0.0;
		sy[i][0] =0.0;
		syold[i][n] =0.0;
		syold[i][0] = 0.0;
	}	
	
	//at left and right wall//
	for(j=0;j<=n;j++)
	{
		sy[0][j] = ((4.0*sy[1][j]) - sy[2][j])/3.0;
		sy[m][j] = ((4.0*sy[m-1][j]) - sy[m-2][j])/3.0;
		syold[0][j] = ((4.0*syold[1][j]) - syold[2][j])/3.0;
		syold[m][j] = ((4.0*sy[m-1][j]) - sy[m-2][j])/3.0;
	}
}

double induced_field_calculation()
{
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			if(j>0&&j<=n1)
			{
				if(i>=0&&i<=m1)
				{
					flag = 0.0;
				}
				else if(i>=m1&&i<=m2)
				{
					flag = 0.0;
				}
				else if(i>=m3&&i<=m4)
				{
					flag = 0.0;
				}
				else if(i>=m5&&i<=m6)
				{
					flag = 0.0;
				}
				else if(i>=m6&&i<=m)
				{
					flag = 0.0;
				}
				else
				{
					flag = 1.0;
				}
				
			}
			else if(j>=n2&&j<n)
			{
				if(i>=m1&&i<=m2)
				{
					flag = 1.0;
				}
				else if(i>=m3&&i<=m4)
				{
					flag = 1.0;
				}
				else if(i>=m5&&i<=m6)
				{
					flag = 1.0;
				}
				else
				{
					flag = 0.0;
				}
			}
			else
			{
				flag = 0.0;
			}
			fakesy[i][j] = (((syold[i+1][j] + sy[i-1][j])/(dx*dx)) + ((syold[i][j+1] + sy[i][j-1])/(dy*dy)) - (kappa*kappa*sinh(syold[i][j])) + (kappap*kappap*flag))/(2.0*((1.0/(dx*dx)) + (1.0/(dy*dy))));
			sy[i][j] = syold[i][j] + (wsy*(fakesy[i][j] - syold[i][j]));
		}		
	}
}

double induced_field_error()
{
	sum = 0.0;
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			sum = sum + ((sy[i][j] - syold[i][j])*(sy[i][j] - syold[i][j]));
		}
	}
	app_error = sqrt(sum/((m-2)*(n-2)));
}

double induced_field_update()
{
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			syold[i][j] = sy[i][j];
		}
	}
}

double induced_field_print()
{
	FILE *fp5;
	fp5 = fopen("induced field.dat","w");
	fprintf(fp5,"zone\tj=%d\ti=%d",m+1,n+1);
	for(i=0;i<=m;i++)
	{
		for(j=0;j<=n;j++)
		{
			fprintf(fp5,"\n%f\t\t%f\t\t%f",x[i][j],y[i][j],sy[i][j]);
		}
	}
	fclose(fp5);
}


//*************************************************************************************************************************************************************************//
																		//CALCULATING FLOW FIELD//
//*************************************************************************************************************************************************************************//
double flow_field()
{
	guess();
	error = 1000.0;
	count = 1;
	while(error>1.0e-13)
	{
		stream_boundary_condition();
		stream_calculation();
		velocity_boundary_condition();
		velocity_calculation();
		vorticity_boundary_condition();
		vorticity_calculation();
		error_calculation();
		printf("\n%d\t%.20f",count,error);	
		count++;
		update();
	}
}

//*************************************************************************************************************************************************************************//
																		//CALCULATING SPECIES TRANSPORT//
//*************************************************************************************************************************************************************************//

double species_transport()
{
	dif = 1/(re*sc);
	app_error = 1000.0;
	count = 1;
	species_initial();
	while(app_error>1.0e-12)
	{
		species_boundary_condition();
		species_calculation();
		species_error();
		printf("\n%d\t\t%.15f",count,app_error);
		count++;
		species_update();
	}
	species_print();
}
double species_initial()
{
	for(i=0;i<=m;i++)
	{
		for(j=0;j<=n;j++)
		{
			sp[i][j] = 0.0;
			spold[i][j] = 0.0;
		}
	}
}
double species_boundary_condition()
{
	//at upper wall and lower wall//
	for(i=0; i<=m; i++)
	{
		sp[i][n] = sp[i][n-1];
		sp[i][0] = sp[i][1];
		spold[i][n] = spold[i][n-1];
		spold[i][0] = spold[i][1];
	}
	
	//at left wall//
	for(j=0; j<(n/2); j++) //lower half inlet
	{
		sp[0][j] = 1.0;
		spold[0][j] = 1.0;
	}
	for(j=(n/2); j<=n; j++) //upper half inlet
	{
		sp[0][j] = 0.0;
		spold[0][j] = 0.0;
	}
	
	//at right wall//
	for(j=0; j<=n; j++)
	{
		sp[m][j] = sp[m-1][j];
		spold[m][j] = spold[m-1][j];
	}
		
}
double species_calculation()
{
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			fakesp[i][j] = (dif*(((spold[i+1][j]+sp[i-1][j])/(dx*dx))+((spold[i][j+1]+sp[i][j-1])/(dy*dy)))
			-(u[i][j]*((spold[i+1][j]-sp[i-1][j])/(2*dx)))-(v[i][j]*((spold[i][j+1]-sp[i][j-1])/(2*dy))))/((2*dif*((1/(dx*dx))+(1/(dy*dy)))));
			sp[i][j] = spold[i][j]+(wsp*(fakesp[i][j] - spold[i][j]));
		}
	}
}		
double species_error()
{
	sum = 0.0;
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			sum = sum + ((sp[i][j] - spold[i][j])*(sp[i][j] - spold[i][j]));
		}
	}
	app_error = sqrt(sum/(m*n));
}
double species_update()
{
	for(i=1;i<m;i++)
	{
		for(j=1;j<n;j++)
		{
			spold[i][j] = sp[i][j];
		}
	}
}
double species_print()
{
	FILE *fp11;
	fp11 = fopen("species transport.dat","w");
	fprintf(fp11,"zone\tj=%d\ti=%d",m+1,n+1);
	for(i=0;i<=m;i++)
	{
		for(j=0;j<=n;j++)
		{
			fprintf(fp11,"\n%f\t\t%f\t\t%f",x[i][j],y[i][j],sp[i][j]);
		}
	}
	fclose(fp11);
}
