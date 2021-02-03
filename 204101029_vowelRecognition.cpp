#include<iostream>
#include<vector>
#include<fstream>
#include<stdlib.h>
#include<cmath>
#include<math.h>
#define N 320    // Frame size
#define p 12     // order of LPC

using namespace std;

long int file_size;  // file size

// this function is used to calculate DC shift on a silence file.(Note : No Dc shift is required in this program).
void nullify_dc_shift(string file)
{
	ifstream input_file;             
	input_file.open(file.c_str());
	long double sum=0,curr,count=0;
	while(!input_file.eof())        // caluclating sum of all the samples in the file.
	{
		input_file>>curr; 
		sum+=curr;
		count++;
	}
	file_size=count;
	cout<<"DC shift : "<<sum/count<<endl;    // taking averge of that sum and that will be the DC shift.
}


// This function is used to normalize the file which is given as argument to the function.
vector<double> normalize(string file_name)
{
	file_size=0;
	vector<double> normalized;  // this array will contain the normalized value of the input file and at last we return this array.
	float curr,max=0;
	float factor;
	ifstream input_file;
	input_file.open(file_name.c_str());
	while(!input_file.eof())  // finding the maximuam absoulute sample in the whole file.
	{
		input_file>>curr;
		if(max<abs(curr))
		max=abs(curr);
		file_size++;
	}
	factor=10000/max;   // then finding the scaling factor w.r.t to 10000
	input_file.close();
	
	input_file.open(file_name.c_str()); 
	while(!input_file.eof())   // scaling all the samples of the file using the sacling factor.
	{
		input_file>>curr;
		curr=curr*factor;
		normalized.push_back(curr);   // then pushing the scaled value in normalized vector.
	}
	input_file.close();
	return normalized;   // return the normalized vector.
}

// This function divides the whole normalized samples in frames of 320 samples.
void framing(vector<double> normalized,vector<vector<double> >&frames)
{
	
	long int i=0,j=0;
	int n=normalized.size()/320;  // n= number of frames in the given file.
		for(i=0;i<n;i++)
		{
			for(int k=0;k<320;k++)
			{
				frames[i][k]=normalized[j++];
			}
		}
}

// this function is used to apply hamming window on each frames which we have generated in previous function.
void hamming(vector<vector<double> >& frames)
{
	long int i,j;
	double w;
	for(i=0;i<frames.size();i++)
	{
		for(j=0;j<320;j++)
		{
			w=0.54-0.46*cos((2*3.141*j)/319);
			frames[i][j]*=w;
		}
	}
}


// this is the function to calculate autocorelation vectors R for the given frames.
void find_R(vector<vector<double> >&frames, vector<vector<double> >&R)
{
	
	
	int i,j,k;
	long double sum=0.0;
	for(i=0;i<frames.size();i++)
	{
		
		for(k=0;k<13;k++)
		{
			sum=0.0;
			for(j=0;j<320-k;j++)
			{
				sum=sum+frames[i][j]*frames[i][j+k];
			}
			
			R[i][k]=sum/320;
		}
	}
}


// This function is for the durbin algorithm
vector<double> durbin(vector<double>r)
{
	vector<double>b(13),an(13);
	double k,E;
	b[0]=1;
	k=-r[1]/r[0];
	b[1]=k;
	E=r[0]*(1-k*k);   // initial value of E i.e. E[0]
	int i,j,s;
	for(i=2;i<13;i++)   // no we have to calculated the remaing values of an
	{
		s=0;
		for(j=1;j<i;j++)
		{
			s+=r[i-j]*b[j];
		}
		s=s+r[i];
		k=-s/E;
		for(j=1;j<i;j++)    //storing the current i-1 values of a in temporary array an so that we don't need any 2-D array for a
		{
			an[j]=b[j]+k*b[i-j];
		}
		
		an[i]=k;
		for(j=1;j<=i;j++)   // updating the array a using the previous value which are stored in temporary array an.
		b[j]=an[j]; 

		E=E*(1-k*k);   //updating the value of E.
	}
	return b;  // returning the vector b which contains the solutions
}


// this function is used to calculate the a_i using the above durbin method
void find_a(vector<vector<double> >R,vector<vector<double> >&a)
{
	vector<double>temp_a(13);  // this a vector will store the temporary solution and then we just update vector a
	long int i,j,k;
	for(i=0;i<R.size();i++)
	{
		
		temp_a=durbin(R[i]);  // applying durbin methode and storing the results in temporary array.
		for(j=0;j<13;j++)  // updating the changes in final array of cepstral coefficients.
		a[i][j]=temp_a[j];
		
	}
}


//this function is used to calculate the cepstral coefficients. using R_i and a_i calculated above.
vector<vector<double> > find_c(vector<vector<double> >R, vector<vector<double> >a)
{
	
	vector<vector<double> >c(R.size(),vector<double>(13));
	long int i,j,m,k;
	double sum=0.0,d;
	for(i=0;i<R.size();i++)   // outer loop : run over each frame in file.
	{
		for(m=1;m<13;m++)     // we need to calculate 12 c_i
		{
			sum=0.0;
			for(k=1;k<=m-1;k++)  // inner loop to calculate the summation part.
			{
				double kf=k,mf=m;
				d=kf/mf;
				sum=sum+d*c[i][k]*a[i][m-k]*(-1);
			}
		
			c[i][m]=sum+a[i][m];
		}
	}
	return c;   // returning the calculated cepstral coefficients.
}


//applying the raised sine window on the cepstral coefficients calculated above.
void raised_sine(vector<vector<double> >&c)
{
	long int i,j;
	double w;
	for(i=0;i<c.size();i++)
	{
		for(j=0;j<13;j++)
		{
			w=(1+6*sin((3.1414*j)/12));  // w[i] calculation
			c[i][j]*=w;
		}
	}
}


// this function is used to make find the average c_i for 10 files of vowel a and storing these c_i of 5 rows and 12 columns in a file.
vector<vector<double> > train_for_a()
{
	double c[50][12];    // array of 50 rows : 5 steady frames of 10 test file = 50 rows, 12 column = 12 cepstral coefficients of each frames
	int k=0;
	ofstream output;
	output.open("final_c_aa.txt");  // an output file where c_i will be stored.
	for(int Q=0;Q<10;Q++)   // loop through all the 10 test files.
	{
		vector<double> normalized;    // normalized value of each 10 test files
		switch(Q)   // normalizing each file one by one in each iteration of for loop.
		{
			case 0: normalized=normalize("a1.txt");break;
			case 1: normalized=normalize("a2.txt");break;
			case 2: normalized=normalize("a3.txt");break;
			case 3: normalized=normalize("a4.txt");break;
			case 4: normalized=normalize("a5.txt");break;
			case 5: normalized=normalize("a6.txt");break;
			case 6: normalized=normalize("a7.txt");break;
			case 7: normalized=normalize("a8.txt");break;
			case 8: normalized=normalize("a9.txt");break;
			case 9: normalized=normalize("a10.txt");break;
		}
		vector<vector<double> >A(normalized.size()/320,vector<double>(13)),C(normalized.size()/320,vector<double>(13));
		vector<vector<double> >frames(normalized.size()/320,vector<double>(320)),R(normalized.size()/320,vector<double>(13));
		framing(normalized,frames);   // make frames of normalized file
		hamming(frames);   // applying hamming window on frames
		find_R(frames,R);  // finding R_i of frames
		find_a(R,A);   // finding a_i for each frames
		C=find_c(R,A);  // finding c_i for each frames
		raised_sine(C);  // applying raised_sin window on c_i
		int marker=0;

		for(long int i=0;i<C.size();i++)
		{
			if(R[i][0]>999999 && marker<5)   // taking 5 steady frames
			{
				for(int j=1;j<13;j++)
				{
					c[k][j-1]=C[i][j];
				}
				k++;
				marker++;
			}
		}	
	}
	
	int i,j;	
	vector<vector<double> >final_c(5,vector<double>(12));
	double sum1,sum2,sum3,sum4,sum5;

	// loop to calculate the average of each 10 test filea and store them in final_c vector of 5 rows 12 cols
	for(j=0;j<12;j++)
	{
			sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;
			sum1=c[0][j]+c[5][j]+c[10][j]+c[15][j]+c[20][j]+c[25][j]+c[30][j]+c[35][j]+c[40][j]+c[45][j];
			sum2=c[0+1][j]+c[5+1][j]+c[10+1][j]+c[15+1][j]+c[20+1][j]+c[25+1][j]+c[30+1][j]+c[35+1][j]+c[40+1][j]+c[45+1][j];
			sum3=c[0+2][j]+c[5+2][j]+c[10+2][j]+c[15+2][j]+c[20+2][j]+c[25+2][j]+c[30+2][j]+c[35+2][j]+c[40+2][j]+c[45+2][j];
			sum4=c[0+3][j]+c[5+3][j]+c[10+3][j]+c[15+3][j]+c[20+3][j]+c[25+3][j]+c[30+3][j]+c[35+3][j]+c[40+3][j]+c[45+3][j];
			sum5=c[0+4][j]+c[5+4][j]+c[10+4][j]+c[15+4][j]+c[20+4][j]+c[25+4][j]+c[30+4][j]+c[35+4][j]+c[40+4][j]+c[45+4][j];
			
			final_c[0][j]=sum1/10;
			final_c[1][j]=sum2/10;
			final_c[2][j]=sum3/10;
			final_c[3][j]=sum4/10;
			final_c[4][j]=sum5/10;
	}
	
	for(i=0;i<5;i++)
	{
		for(j=0;j<12;j++)
		{
			output<<final_c[i][j]<<" ";
		}
		output<<endl;
	}
	output.close();
	return final_c;
}



// this function will work same as train_for_a() but it will train the values of vowel e.
vector<vector<double> > train_for_e()
{
	double c[50][12];
	int k=0;
	ofstream output;
	output.open("final_c_ee.txt");
	for(int Q=0;Q<10;Q++)
	{
		vector<double> normalized;
		switch(Q)
		{
			case 0: normalized=normalize("e1.txt");break;
			case 1: normalized=normalize("e2.txt");break;
			case 2: normalized=normalize("e3.txt");break;
			case 3: normalized=normalize("e4.txt");break;
			case 4: normalized=normalize("e5.txt");break;
			case 5: normalized=normalize("e6.txt");break;
			case 6: normalized=normalize("e7.txt");break;
			case 7: normalized=normalize("e8.txt");break;
			case 8: normalized=normalize("e9.txt");break;
			case 9: normalized=normalize("e10.txt");break;
		}
		vector<vector<double> >A(normalized.size()/320,vector<double>(13)),C(normalized.size()/320,vector<double>(13));
		vector<vector<double> >frames(normalized.size()/320,vector<double>(320)),R(normalized.size()/320,vector<double>(13));
		framing(normalized,frames);
		hamming(frames);
		find_R(frames,R);
		find_a(R,A);
		C=find_c(R,A);
		raised_sine(C);
		int marker=0;
		for(long int i=0;i<C.size();i++)
		{
			if(R[i][0]>999999 && marker<5)
			{
				for(int j=1;j<13;j++)
				{
					c[k][j-1]=C[i][j];
				}
				k++;
				marker++;
			}
		}	
	}
	
	int i,j;	
	vector<vector<double> >final_c(5,vector<double>(12));
	double sum1,sum2,sum3,sum4,sum5;
	for(j=0;j<12;j++)
	{
			sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;
			sum1=c[0][j]+c[5][j]+c[10][j]+c[15][j]+c[20][j]+c[25][j]+c[30][j]+c[35][j]+c[40][j]+c[45][j];
			sum2=c[0+1][j]+c[5+1][j]+c[10+1][j]+c[15+1][j]+c[20+1][j]+c[25+1][j]+c[30+1][j]+c[35+1][j]+c[40+1][j]+c[45+1][j];
			sum3=c[0+2][j]+c[5+2][j]+c[10+2][j]+c[15+2][j]+c[20+2][j]+c[25+2][j]+c[30+2][j]+c[35+2][j]+c[40+2][j]+c[45+2][j];
			sum4=c[0+3][j]+c[5+3][j]+c[10+3][j]+c[15+3][j]+c[20+3][j]+c[25+3][j]+c[30+3][j]+c[35+3][j]+c[40+3][j]+c[45+3][j];
			sum5=c[0+4][j]+c[5+4][j]+c[10+4][j]+c[15+4][j]+c[20+4][j]+c[25+4][j]+c[30+4][j]+c[35+4][j]+c[40+4][j]+c[45+4][j];
			
			final_c[0][j]=sum1/10;
			final_c[1][j]=sum2/10;
			final_c[2][j]=sum3/10;
			final_c[3][j]=sum4/10;
			final_c[4][j]=sum5/10;
	}
	
	for(i=0;i<5;i++)
	{
		for(j=0;j<12;j++)
		{
			output<<final_c[i][j]<<" ";
		}
		output<<endl;
	}
	output.close();
	return final_c;
}



// same as above:  train_for_e(), train_for_a()
vector<vector<double> > train_for_u()
{
	double c[50][12];
	int k=0;
	ofstream output;
	output.open("final_c_u.txt");
	for(int Q=0;Q<10;Q++)
	{
		vector<double> normalized;
		switch(Q)
		{
			case 0: normalized=normalize("u1.txt");break;
			case 1: normalized=normalize("u2.txt");break;
			case 2: normalized=normalize("u3.txt");break;
			case 3: normalized=normalize("u4.txt");break;
			case 4: normalized=normalize("u5.txt");break;
			case 5: normalized=normalize("u6.txt");break;
			case 6: normalized=normalize("u7.txt");break;
			case 7: normalized=normalize("u8.txt");break;
			case 8: normalized=normalize("u9.txt");break;
			case 9: normalized=normalize("u10.txt");break;
		}
		vector<vector<double> >A(normalized.size()/320,vector<double>(13)),C(normalized.size()/320,vector<double>(13));
		vector<vector<double> >frames(normalized.size()/320,vector<double>(320)),R(normalized.size()/320,vector<double>(13));
		framing(normalized,frames);
		hamming(frames);
		find_R(frames,R);
		find_a(R,A);
		C=find_c(R,A);
		raised_sine(C);
		int marker=0;
		for(long int i=0;i<C.size();i++)
		{
			if(R[i][0]>999999 && marker<5)
			{
				for(int j=1;j<13;j++)
				{
					c[k][j-1]=C[i][j];
				}
				k++;
				marker++;
			}
		}	
	}
	
	int i,j;	
	vector<vector<double> >final_c(5,vector<double>(12));
	double sum1,sum2,sum3,sum4,sum5;
	for(j=0;j<12;j++)
	{
			sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;
			sum1=c[0][j]+c[5][j]+c[10][j]+c[15][j]+c[20][j]+c[25][j]+c[30][j]+c[35][j]+c[40][j]+c[45][j];
			sum2=c[0+1][j]+c[5+1][j]+c[10+1][j]+c[15+1][j]+c[20+1][j]+c[25+1][j]+c[30+1][j]+c[35+1][j]+c[40+1][j]+c[45+1][j];
			sum3=c[0+2][j]+c[5+2][j]+c[10+2][j]+c[15+2][j]+c[20+2][j]+c[25+2][j]+c[30+2][j]+c[35+2][j]+c[40+2][j]+c[45+2][j];
			sum4=c[0+3][j]+c[5+3][j]+c[10+3][j]+c[15+3][j]+c[20+3][j]+c[25+3][j]+c[30+3][j]+c[35+3][j]+c[40+3][j]+c[45+3][j];
			sum5=c[0+4][j]+c[5+4][j]+c[10+4][j]+c[15+4][j]+c[20+4][j]+c[25+4][j]+c[30+4][j]+c[35+4][j]+c[40+4][j]+c[45+4][j];
			
			final_c[0][j]=sum1/10;
			final_c[1][j]=sum2/10;
			final_c[2][j]=sum3/10;
			final_c[3][j]=sum4/10;
			final_c[4][j]=sum5/10;
	}
	
	for(i=0;i<5;i++)
	{
		for(j=0;j<12;j++)
		{
			output<<final_c[i][j]<<" ";
		}
		output<<endl;
	}
	output.close();
	return final_c;
}




vector<vector<double> > train_for_i()
{
	double c[50][12];
	int k=0;
	ofstream output;
	output.open("final_c_i.txt");
	for(int Q=0;Q<10;Q++)
	{
		vector<double> normalized;
		switch(Q)
		{
			case 0: normalized=normalize("i1.txt");break;
			case 1: normalized=normalize("i2.txt");break;
			case 2: normalized=normalize("i3.txt");break;
			case 3: normalized=normalize("i4.txt");break;
			case 4: normalized=normalize("i5.txt");break;
			case 5: normalized=normalize("i6.txt");break;
			case 6: normalized=normalize("i7.txt");break;
			case 7: normalized=normalize("i8.txt");break;
			case 8: normalized=normalize("i9.txt");break;
			case 9: normalized=normalize("i10.txt");break;
		}
		vector<vector<double> >A(normalized.size()/320,vector<double>(13)),C(normalized.size()/320,vector<double>(13));
		vector<vector<double> >frames(normalized.size()/320,vector<double>(320)),R(normalized.size()/320,vector<double>(13));
		framing(normalized,frames);
		hamming(frames);
		find_R(frames,R);
		find_a(R,A);
		C=find_c(R,A);
		raised_sine(C);
		int marker=0;
		for(long int i=0;i<C.size();i++)
		{
			if(R[i][0]>999999 && marker<5)
			{
				for(int j=1;j<13;j++)
				{
					c[k][j-1]=C[i][j];
				}
				k++;
				marker++;
			}
		}	
	}
	
	int i,j;	
	vector<vector<double> >final_c(5,vector<double>(12));
	double sum1,sum2,sum3,sum4,sum5;
	for(j=0;j<12;j++)
	{
			sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;
			sum1=c[0][j]+c[5][j]+c[10][j]+c[15][j]+c[20][j]+c[25][j]+c[30][j]+c[35][j]+c[40][j]+c[45][j];
			sum2=c[0+1][j]+c[5+1][j]+c[10+1][j]+c[15+1][j]+c[20+1][j]+c[25+1][j]+c[30+1][j]+c[35+1][j]+c[40+1][j]+c[45+1][j];
			sum3=c[0+2][j]+c[5+2][j]+c[10+2][j]+c[15+2][j]+c[20+2][j]+c[25+2][j]+c[30+2][j]+c[35+2][j]+c[40+2][j]+c[45+2][j];
			sum4=c[0+3][j]+c[5+3][j]+c[10+3][j]+c[15+3][j]+c[20+3][j]+c[25+3][j]+c[30+3][j]+c[35+3][j]+c[40+3][j]+c[45+3][j];
			sum5=c[0+4][j]+c[5+4][j]+c[10+4][j]+c[15+4][j]+c[20+4][j]+c[25+4][j]+c[30+4][j]+c[35+4][j]+c[40+4][j]+c[45+4][j];
			
			final_c[0][j]=sum1/10;
			final_c[1][j]=sum2/10;
			final_c[2][j]=sum3/10;
			final_c[3][j]=sum4/10;
			final_c[4][j]=sum5/10;
	}
	
	for(i=0;i<5;i++)
	{
		for(j=0;j<12;j++)
		{
			output<<final_c[i][j]<<" ";
		}
		output<<endl;
	}
	output.close();
	return final_c;
}



vector<vector<double> > train_for_o()
{
	double c[50][12];
	int k=0;
	ofstream output;
	output.open("final_c_o.txt");
	for(int Q=0;Q<10;Q++)
	{
		vector<double> normalized;
		switch(Q)
		{
			case 0: normalized=normalize("o1.txt");break;
			case 1: normalized=normalize("o2.txt");break;
			case 2: normalized=normalize("o3.txt");break;
			case 3: normalized=normalize("o4.txt");break;
			case 4: normalized=normalize("o5.txt");break;
			case 5: normalized=normalize("o6.txt");break;
			case 6: normalized=normalize("o7.txt");break;
			case 7: normalized=normalize("o8.txt");break;
			case 8: normalized=normalize("o9.txt");break;
			case 9: normalized=normalize("o10.txt");break;
		}
		vector<vector<double> >A(normalized.size()/320,vector<double>(13)),C(normalized.size()/320,vector<double>(13));
		vector<vector<double> >frames(normalized.size()/320,vector<double>(320)),R(normalized.size()/320,vector<double>(13));
		framing(normalized,frames);
		hamming(frames);
		find_R(frames,R);
		find_a(R,A);
		C=find_c(R,A);
		raised_sine(C);
		int marker=0;
		for(long int i=0;i<C.size();i++)
		{
			if(R[i][0]>999999 && marker<5)
			{
				for(int j=1;j<13;j++)
				{
					c[k][j-1]=C[i][j];
				}
				k++;
				marker++;
			}
		}	
	}
	
	int i,j;	
	vector<vector<double> >final_c(5,vector<double>(12));
	double sum1,sum2,sum3,sum4,sum5;
	for(j=0;j<12;j++)
	{
			sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;
			sum1=c[0][j]+c[5][j]+c[10][j]+c[15][j]+c[20][j]+c[25][j]+c[30][j]+c[35][j]+c[40][j]+c[45][j];
			sum2=c[0+1][j]+c[5+1][j]+c[10+1][j]+c[15+1][j]+c[20+1][j]+c[25+1][j]+c[30+1][j]+c[35+1][j]+c[40+1][j]+c[45+1][j];
			sum3=c[0+2][j]+c[5+2][j]+c[10+2][j]+c[15+2][j]+c[20+2][j]+c[25+2][j]+c[30+2][j]+c[35+2][j]+c[40+2][j]+c[45+2][j];
			sum4=c[0+3][j]+c[5+3][j]+c[10+3][j]+c[15+3][j]+c[20+3][j]+c[25+3][j]+c[30+3][j]+c[35+3][j]+c[40+3][j]+c[45+3][j];
			sum5=c[0+4][j]+c[5+4][j]+c[10+4][j]+c[15+4][j]+c[20+4][j]+c[25+4][j]+c[30+4][j]+c[35+4][j]+c[40+4][j]+c[45+4][j];
			
			final_c[0][j]=sum1/10;
			final_c[1][j]=sum2/10;
			final_c[2][j]=sum3/10;
			final_c[3][j]=sum4/10;
			final_c[4][j]=sum5/10;
	}
	
	for(i=0;i<5;i++)
	{
		for(j=0;j<12;j++)
		{
			output<<final_c[i][j]<<" ";
		}
		output<<endl;
	}
	output.close();
	return final_c;
}




// this function is used to make a test file ready in the sense that apply normalization, framing, etc 
vector<vector<double> > get_test(string test_file)
{
	    vector<double>normalized;
		normalized=normalize(test_file.c_str());  // normalizing the test file
		vector<vector<double> >test(5,vector<double>(12));
		vector<vector<double> >A(normalized.size()/320,vector<double>(13)),C(normalized.size()/320,vector<double>(13));
		vector<vector<double> >frames(normalized.size()/320,vector<double>(320)),R(normalized.size()/320,vector<double>(13));
		framing(normalized,frames);     // framing the test file
		hamming(frames);  // applying the hamming window
		find_R(frames,R);  // finding the R_i
		find_a(R,A);    // finding the A_i
		C=find_c(R,A);  // finding the C_i
		raised_sine(C);  // raised sine window
		int marker=0,k=0; 
		//finding 5 steady frames from test file
		for(long int i=0;i<C.size();i++)
		{
			if(R[i][0]>999999 && marker<5)
			{
				for(int j=1;j<13;j++)
				{
					test[k][j-1]=C[i][j];
				}
				k++;
				marker++;
			}
		}	
	return test;
}

// reference repository of each vowel
vector<vector<double> >reference_a(5,vector<double>(12));
vector<vector<double> >reference_e(5,vector<double>(12));
vector<vector<double> >reference_i(5,vector<double>(12));
vector<vector<double> >reference_o(5,vector<double>(12));
vector<vector<double> >reference_u(5,vector<double>(12));

//this function will predict the vowel by calculating the tokhura's distance of test file from each refrence file
char makeprediction(vector<vector<double> >test)
{
	long double curr=0.0,min;
	double w[12]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};  //tokhura's weights given by sir
	
	char chr='a';  // initially assume that predicted charater is 'a'
	int i,j;
	for(i=0;i<5;i++)
	{
		for(j=0;j<12;j++)
		{
			curr=curr+w[j]*pow((reference_a[i][j]-test[i][j]),2);   // calculating the tokhura's distance of test file from reference file of a
		}
	} 
	min=curr;
	curr=0.0;
	for(i=0;i<5;i++)
	{
		for(j=0;j<12;j++)
		{
			curr=curr+w[j]*pow((reference_e[i][j]-test[i][j]),2);// calculating the tokhura's distance of test file from reference file of e.
		}
	}
	if(min>curr) // if tokhura's distace of e is less than that of a.
	{
		min=curr;  // assume new min as tokhura's distance of e
		chr='e';  // predicted char is e
	}
	curr=0.0;
	for(i=0;i<5;i++)
	{
		for(j=0;j<12;j++)
		{
			curr=curr+w[j]*pow((reference_i[i][j]-test[i][j]),2);  // calculating the tokhura's distance of test file from reference file of i
		}
	}
	
	if(min>curr)  // if tokhura's distance of i is less than the current minimum
	{ 
		min=curr;  // new minimum is current minimum
		chr='i';  // new prediction is i
	}
	curr=0.0;
	for(i=0;i<5;i++)
	{
		for(j=0;j<12;j++)
		{
			curr=curr+w[j]*pow((reference_o[i][j]-test[i][j]),2);
		}
	}
	
	if(min>curr)
	{
		min=curr;
		chr='o';
	}
	curr=0.0;
	for(i=0;i<5;i++)
	{
		for(j=0;j<12;j++)
		{
			curr=curr+w[j]*pow((reference_u[i][j]-test[i][j]),2);
		}
	}
	if(min>curr)
	{
		min=curr;
		chr='u';
	}
	curr=0.0;
	
	return chr;  //returning the predicted character
	
}


//main drivr function
int main()
{
	cout<<"-----------------------------------> CS566:Speech Processing, Assignemnt 3 <------------------------------------\n";
	cout<<"Submitted To:Dr. Pradip Kumar Das, TA's:Mayank sharma,Komal bharti,Vandana mishra,Sakshi sharma,Prabattina bhath\n";
    cout<<"-----------------------------------> Submitted By : Himanshu Patel (204101029)<----------------------------------\n";
    cout<<"*****************************************************************************************************************\n\n";
    cout<<"Training....\n";
	ofstream output;
    char chr;
    output.open("prediction.txt");
	reference_a=train_for_a();  // generating the reference file for a
	reference_e=train_for_e();  // generating the reference file for e
	reference_u=train_for_u();  // generating the reference file for i
	reference_i=train_for_i();  // generating the reference file for o
	reference_o=train_for_o();  // generating the reference file for u
	cout<<"Training completed!!\n";
	vector<vector<double> >test(5,vector<double> (12));
	int n=1;
	int i,j;
	while(true)
	{
		cout<<"Enter 0 for live recording test, 1 for mannual test, -1 to exit\n";
		cin>>n;
		if(n==-1)
		exit(0);
		else if(n==1)
		{
			cout<<"Enter the file number to be tested enter 1-10 for a, 11 to 20 for e, 21 to 30 for i,31 to 40 for o, 41 to 50 for u\n";
			cin>>n;
		}
		switch(n)
		{
			case 0: system("Recording_Module.exe 3 live_recording.wav live_recording.txt"); test=get_test("live_recording.txt");break;
			case 1:test=get_test("a11.txt");break;
			case 2:test=get_test("a12.txt");break;
			case 3:test=get_test("a13.txt");break;
			case 4:test=get_test("a14.txt");break;
			case 5:test=get_test("a15.txt");break;
			case 6:test=get_test("a16.txt");break;
			case 7:test=get_test("a17.txt");break;
			case 8:test=get_test("a18.txt");break;
			case 9:test=get_test("a19.txt");break;
			case 10:test=get_test("a20.txt");break;
			case 11:test=get_test("e11.txt");break;
			case 12:test=get_test("e12.txt");break;
			case 13:test=get_test("e13.txt");break;
			case 14:test=get_test("e14.txt");break;
			case 15:test=get_test("e15.txt");break;
			case 16:test=get_test("e16.txt");break;
			case 17:test=get_test("e17.txt");break;
			case 18:test=get_test("e18.txt");break;
			case 19:test=get_test("e19.txt");break;
			case 20:test=get_test("e20.txt");break;
			case 21:test=get_test("i11.txt");break;
			case 22:test=get_test("i12.txt");break;
			case 23:test=get_test("i13.txt");break;
			case 24:test=get_test("i14.txt");break;
			case 25:test=get_test("i15.txt");break;
			case 26:test=get_test("i16.txt");break;
			case 27:test=get_test("i17.txt");break;
			case 28:test=get_test("i18.txt");break;
			case 29:test=get_test("i19.txt");break;
			case 30:test=get_test("i10.txt");break;
			case 31:test=get_test("o11.txt");break;
			case 32:test=get_test("o12.txt");break;
			case 33:test=get_test("o13.txt");break;
			case 34:test=get_test("o14.txt");break;
			case 35:test=get_test("o15.txt");break;
			case 36:test=get_test("o16.txt");break;
			case 37:test=get_test("o17.txt");break;
			case 38:test=get_test("o18.txt");break;
			case 39:test=get_test("o19.txt");break;
			case 40:test=get_test("o20.txt");break;
			case 41:test=get_test("u11.txt");break;
			case 42:test=get_test("u12.txt");break;
			case 43:test=get_test("u13.txt");break;
			case 44:test=get_test("u14.txt");break;
			case 45:test=get_test("u15.txt");break;
			case 46:test=get_test("u16.txt");break;
			case 47:test=get_test("u17.txt");break;
			case 48:test=get_test("u18.txt");break;
			case 49:test=get_test("u19.txt");break;
			case 50:test=get_test("u20.txt");break;
			case 51: exit(0);
			default : cout<<"Invalid\n"; break;

		}
		
		chr=makeprediction(test);
		cout<<"Predicted character is : "<<chr;
		cout<<endl;
		cout<<"*******************************************************************************\n\n";
	}
	output.close();
	return 0;
}

