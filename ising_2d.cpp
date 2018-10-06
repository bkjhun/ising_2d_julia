#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MOD(a,b) ((((a)%(b))+(b))%(b))

int main(int argc, char* argv[])
{
	int Length = atoi(argv[1]);
	float Temperature;
	int max_index_Temperature;

	int num_sample = 10000, num_sample_burn = 100;
	int period_sample;
	
	int E0; 
	double magnetization_sample, magnetization, magnetization2;
	double energy_sample, energy, energy2;
	int x,y;
	char FilenameOut[100];
	float vec_Temperature[100];
	double list_exp[3];

	FILE* FileOut;

	FileOut = fopen("./data/list_Temperature.txt", "r");
	for(int index_Temperature=0; index_Temperature < 100; index_Temperature++)
	{
		if(feof(FileOut))
		{
			max_index_Temperature = index_Temperature;
			break;
		}
		fscanf(FileOut, "%f\n", &vec_Temperature[index_Temperature]);
	}
	fclose(FileOut);
	
	sprintf(FilenameOut, "./data/magnetization_L%d.txt", Length);
	FileOut = fopen(FilenameOut, "w");
	
	char *SpinConfig = (char*)malloc(sizeof(char)*Length*Length);

	for(int num=0; num<Length*Length; num++)
	{
		SpinConfig[num] = 2*(rand()%2) -1;
		//SpinConfig[num] = 1;
	}

	for(int index_Temperature=0; index_Temperature < max_index_Temperature; index_Temperature++)
	{
		Temperature = vec_Temperature[index_Temperature];
		
		if(Temperature > 2.1 && Temperature < 2.5)
			period_sample = 30 *Length*Length /16/16;
		else
			period_sample = 30;
		
		list_exp[0] = exp(-8.0/Temperature);
		list_exp[1] = exp(-4.0/Temperature);
		
		for(int index_sample=0; index_sample<num_sample_burn*Length*Length; index_sample++)
		{
			x = rand()%Length;
			y = rand()%Length;
			
			E0 = SpinConfig[x*Length+y] * (SpinConfig[Length*MOD(x-1,Length) + y] + SpinConfig[Length*MOD(x+1,Length) + y] + SpinConfig[Length*x + MOD(y-1,Length)] + SpinConfig[Length*x + MOD(y+1,Length)]);
			
			if(E0 >= 0)
				SpinConfig[x*Length+y] = -SpinConfig[x*Length+y];
			else if(rand()/(float)RAND_MAX < list_exp[E0/2+2])
				SpinConfig[x*Length+y] = -SpinConfig[x*Length+y];
		}

		magnetization = 0;
		magnetization2 = 0;
		energy = 0;
		energy2 = 0;
		
		for(int index_sample=0; index_sample<num_sample; index_sample++)
		{
			for(int iteration=0; iteration<period_sample; iteration++)
			{
				for(int x=0; x<Length; x++)
				{
					for(int y=0; y<Length; y++)
					{
						E0 = -SpinConfig[x*Length+y] * (SpinConfig[Length*MOD(x-1,Length) + y] + SpinConfig[Length*MOD(x+1,Length) + y] + SpinConfig[Length*x + MOD(y-1,Length)] + SpinConfig[Length*x + MOD(y+1,Length)]);
						
						if(E0 >= 0)
							SpinConfig[x*Length+y] = -SpinConfig[x*Length+y];
						else if(rand()/(float)RAND_MAX < list_exp[E0/2+2])
							SpinConfig[x*Length+y] = -SpinConfig[x*Length+y];						
					}
				}
			}

			magnetization_sample = 0;
			energy_sample = 0;
			for(x=0; x<Length; x++)
			{
				for(y=0; y<Length; y++)
				{
					magnetization_sample += SpinConfig[x*Length+y];
					energy_sample += -SpinConfig[x*Length+y] * (SpinConfig[Length*MOD(x-1,Length) + y] + SpinConfig[Length*MOD(x+1,Length) + y] + SpinConfig[Length*x + MOD(y-1,Length)] + SpinConfig[Length*x + MOD(y+1,Length)]);
				}
			}
			magnetization += fabs(magnetization_sample) /(float)num_sample /Length/Length;
			magnetization2 += magnetization_sample * magnetization_sample /(float)num_sample /Length/Length/Length/Length;
			energy += energy_sample /(float)num_sample /Length/Length;
			energy2 += energy_sample * energy_sample /(float)num_sample /Length/Length/Length/Length;
		}
		fprintf(FileOut, "%lf\t%lf\t%lf\t%lf\t%lf\n", Temperature, magnetization, magnetization2-magnetization*magnetization, energy, energy2-energy*energy);
		fflush(FileOut);

	}
}
