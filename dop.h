/***************************************************************************************************************************************************\
*								 Class: DOP					 				   																		*
* Tinos & Yang (2014). "Analysis of fitness landscape modifications in evolutionary dynamic optimization", Information Sciences, 282.				*
* DOP types:																																		*
* 	DOP Type 1: DOP with permutation																												*
* 		DOP Type 1.1 (DOPs with permutation of the candidate solutions defined by candidate solution exchanges of the XOR type)						*
* 		DOP Type 1.2 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a permutation matrix)	*
* 		DOP Type 1.3 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a set of templates)		*
* 	DOP Type 2: Single time-dependent DOP obtained by duplication of decision variables																*
* 		DOP Type 2.1 (DOPs obtained by copying elements of the decision variables according to a linear transformation)								*
* 		DOP Type 2.2 (DOPs obtained by copying decision variables according to a set of templates)													*
* 	DOP Type 3: Single time-dependent DOPs obtained by adding fitness terms according to a set of templates											*			
\***************************************************************************************************************************************************/

/***************************************************************************************************************************************************\
* Copyright (C) 2016  Renato Tinos & Shengxiang Yang
 * 
 * Reference:
 *		1) Tinos, R. & Yang, S. (2014). "Analysis of fitness landscape modifications in evolutionary dynamic optimization", Information Sciences, 282: 214-236.
 *
 * dop.h is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * dop.h  is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
\***************************************************************************************************************************************************/

#include <cmath>
#include <cstdlib>
#define N_max_s 3

class dop{
	private:	
		int l;															// size of the chromosome
		int change_type;												// change type: 0: no change 1:type 1.1 , 2:type 1.2 , 3:type 1.3 , 4:type 2.1 , 5:type 2.2 , 6:type 3
		int flag_t13;													// flag indicating that type 1.3 is in the first cycle
		int ns;															// number of templates (schemata)
		int os;															// order of the templates (schemata)
		int *m_t11, *b_t12, *l_t21, *m_t22, *s_t22, df_t30[N_max_s];	// control parameters
		int **M_t13, **S_t13, **S_t30;									// control parameters
		void resetDOP(void);											// reset DOP generator
	public:		
		dop(int l_p);		
	    ~dop(void);
		double transform( int *x ,  int *x_new);						// fitness transformation
		void change( int change_type_p, double rho, double f_range ); 	// change in the DOP
}; 



/******************************************************************************\
*							Reset DOP generator								   *
\******************************************************************************/
void dop::resetDOP(void){
	int i, j;
		
	// Initial values for the control parameters
	for (i=0;i<l;i++)
		m_t11[i]=0;
	for (i=0;i<l;i++)
		b_t12[i]=i;
	for (i=0;i<l;i++)
		for (j=0;j<N_max_s;j++)	
	  		S_t13[i][j]=0;
	flag_t13=0;
	for (i=0;i<l;i++)
		l_t21[i]=i;
	for (i=0;i<l;i++)
		s_t22[i]=0;

}


/******************************************************************************\
*								Constructor													   *
\******************************************************************************/
dop::dop(int l_p){
	
	// Parameters
	l=l_p;				// dimension (solution vector lenght)
	
	// Memory allocation 
	m_t11=aloc_vectori(l);
	b_t12=aloc_vectori(l);
	l_t21=aloc_vectori(l);
	m_t22=aloc_vectori(l);
	s_t22=aloc_vectori(l);
	M_t13=aloc_matrixi(l,N_max_s);
	S_t13=aloc_matrixi(l,N_max_s);
	S_t30=aloc_matrixi(l,N_max_s);
	
	change_type=0;		// stationary landscape
	resetDOP();
		
}


/******************************************************************************\
*								 Destructor													   *
\******************************************************************************/
dop::~dop(void){

	delete [] m_t11;
	delete [] b_t12;
	delete [] l_t21;
	delete [] m_t22;
	delete [] s_t22;
	desaloc_matrixi(M_t13,l);
	desaloc_matrixi(S_t13,l);
  	desaloc_matrixi(S_t30,l);
	               
}  


/******************************************************************************\
*								 Fitness transformation			 			   *
\******************************************************************************/
double dop::transform( int *x ,  int *x_new){
	int i, k, schema_flag, *m_temp;
	double df;
	
	if (change_type==0){
		// No change
		for (i=0;i<l;i++)
			x_new[i]=x[i];
		df=0.0;		
	} 
	else if (change_type==1){
		// DOP Type 1.1 (DOPs with permutation of the candidate solutions defined by candidate solution exchanges of the XOR type)	
		XOR(x,m_t11,x_new,l);						// Eq. 35
		df=0.0;
	} 
	else if (change_type==2){
		// DOP Type 1.2 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a permutation matrix)
		// Eq. 37: observation - we don´t use a matrix B_t12 because of it costs O(l^2)
		for (i=0;i<l;i++)
			x_new[i]=x[b_t12[i]];
		df=0.0;
	} 
	else if (change_type==3){
		// DOP Type 1.3 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a set of templates)
		for (i=0;i<l;i++)
			x_new[i]=x[i];
		// Comparing x with the ns templates and changing the respective bits
		for (k=0;k<ns;k++){	
            schema_flag=1 ;           // 0 if schema m is not present in x and 1 otherwise
            i=0;
            while (schema_flag>0 && i<l){
                if (S_t13[i][k]>=0)
                    if (S_t13[i][k] != x[i] )
                        schema_flag=0;
                i++;
        	}
        	// Eq. 39
            if (schema_flag>0){
            	m_temp=aloc_vectori(l);
            	for (i=0;i<l;i++)
            		m_temp[i]=M_t13[i][k];
				XOR(x,m_temp,x_new,l);  
				delete [] m_temp;                
       		}	
		} 
		df=0.0;
	}
	else if (change_type==4){
		// DOP Type 2.1 (DOPs obtained by copying elements of the decision variables according to a linear transformation)	
		// Eq. 43: observation - we don´t use a matrix L_t21 because of it costs O(l^2)
		for (i=0;i<l;i++)
			x_new[i]=x[l_t21[i]];
		df=0.0;
	} 
	else if (change_type==5){
		// DOP Type 2.2 (DOPs obtained by copying decision variables according to a set of templates)															
		// Comparing x with the template and changing the respective bits
	    schema_flag=1 ;           // 0 if schema m is not present in x and 1 otherwise
        i=0;
        while (schema_flag>0 && i<l){
        	if (s_t22[i]>=0)
            	if (s_t22[i] != x[i] )
                	schema_flag=0;
            i++;
        }
        // Eq. 45
        if (schema_flag>0){        
			for (i=0;i<l;i++)
				x_new[i]=m_t22[i];
		}
		else {		
			for (i=0;i<l;i++)
				x_new[i]=x[i];	
		}			
		df=0.0;
	} 
	else if (change_type==6){
		// DOP Type 3: Single time-dependent DOPs obtained by adding fitness terms according to a set of templates	
		for (i=0;i<l;i++)
			x_new[i]=x[i];
		df=0.0;
		// Comparing x with the ns templates and changing the respective bits
		for (k=0;k<ns;k++){	
            schema_flag=1 ;           // 0 if schema m is not present in x and 1 otherwise
            i=0;
            while (schema_flag>0 && i<l){
                if (S_t30[i][k]>=0)
                    if (S_t30[i][k] != x[i] )
                        schema_flag=0;
                i++;
        	}
        	// Eq. 47, 48, 49
            if (schema_flag>0)
            	df=df+df_t30[k];
        }            	
	} 
	
	return df;
}


/******************************************************************************\
*								Change in the DOP							   *
\******************************************************************************/
void dop::change( int change_type_p, double rho, double f_range ){

	int i, k, l_minus_order2, temp, line_1, line_2, *perm_v, *temp_v, *r;
	double df_range;
	
	change_type=change_type_p;

	// random permutation of a vector of integers
	perm_v=aloc_vectori(l);			
	temp_v=aloc_vectori(l);			
	for (i=0;i<l;i++)
		temp_v[i]=i;

 	// number of schemata 
    if (rho <=0.625)   
     	ns=1; 
    else if (rho >0.625 && rho<=0.8125 )       
    	ns=2;
    else     
    	ns=3; 
    // order of the schemata 
    if (rho >= 0.375)        
    	os=1;
    else if (rho <0.375 && rho>=0.1875 )       
    	os=2;
    else       
    	os=3;            
    l_minus_order2 = (l-os)/2;          // the order of a schema is the number of fixed positions in the schema 

			
	if (change_type==0){
		// No change
		resetDOP();			// reset DOP generator
	} 
	else if (change_type==1){
		// DOP Type 1.1 (DOPs with permutation of the candidate solutions defined by candidate solution exchanges of the XOR type)	
		rand_perm(temp_v,perm_v,l);
		r=aloc_vectori(l);
		i=0;
		while (i<l){		
			if(i<rho*l)
				r[perm_v[i]]=1;
			else
				r[perm_v[i]]=0;
			i++;
		}
		for (i=0;i<l;i++)
			temp_v[i]=m_t11[i];
		XOR(temp_v,r,m_t11,l);		// Eq. 36		
		delete [] r;	
	} 
	else if (change_type==2){
		// DOP Type 1.2 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a permutation matrix)
		// permuting 2 lines of matrix B_t12
		// Eq. 38
		for (i=0;i<rho*l;i++){		
			line_1=random_int(0,l-1);
			line_2=random_int(0,l-1);
			while (line_2==line_1)
				line_2=random_int(0,l-1);							
			temp=b_t12[line_1];			
			b_t12[line_1]=b_t12[line_2];
			b_t12[line_2]=temp;			
		}
	} 
	else if (change_type==3){
		// DOP Type 1.3 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a set of templates)
		// Eq. 40
		if (flag_t13==0){
			for (k=0;k<ns;k++){			
				rand_perm(temp_v,perm_v,l);
				i=0;
				while (i<l){		
					if(i<os){					
						S_t13[perm_v[i]][k]=random_int(0,1);
						M_t13[perm_v[i]][k]=0;
					}
					else {					
						S_t13[perm_v[i]][k]=-1;				// -1 indicates don´t care bits
						if(i>l_minus_order2)
							M_t13[perm_v[i]][k]=1;
						else
							M_t13[perm_v[i]][k]=0;
					}
					i++;
				}
			}
		}				
		else{
			for (i=0;i<os;i++){		
				line_1=random_int(0,l-1);
				line_2=random_int(0,l-1);
				while (line_2==line_1)
					line_2=random_int(0,l-1);			
				for (k=0;k<ns;k++){
					// permuting S_t13 and also M_t13		
					temp=S_t13[line_1][k];			
					S_t13[line_1][k]=S_t13[line_2][k];
					S_t13[line_2][k]=temp;
					temp=M_t13[line_1][k];			
					M_t13[line_1][k]=M_t13[line_2][k];
					M_t13[line_2][k]=temp;
				}
			}
		}
		flag_t13=1;
	}
	else if (change_type==4){
		// DOP Type 2.1 (DOPs obtained by copying elements of the decision variables according to a linear transformation)	
		// Eq. 44
		for (i=0;i<l;i++)
			l_t21[i]=i;
		temp=(int) floor(rho*l/2.0);
		for (i=0;i<temp;i++){		
			line_1=random_int(0,l-1);
			line_2=random_int(0,l-1);
			while (line_2==line_1)
				line_2=random_int(0,l-1);	
			l_t21[line_1]=l_t21[line_2];
		}
	} 
	else if (change_type==5){
		// DOP Type 2.2 (DOPs obtained by copying decision variables according to a set of templates)
		// Eq. 46															
		os=l-(int) floor (rho*l/2.0);
		rand_perm(temp_v,perm_v,l);
		i=0;
		while (i<l){		
			if(i<os){					
				s_t22[perm_v[i]]=random_int(0,1);
				m_t22[perm_v[i]]=s_t22[perm_v[i]];
			}
			else {					
				s_t22[perm_v[i]]=-1;				// -1 indicates don´t care bits
				m_t22[perm_v[i]]=random_int(0,1);
					
			}
			i++;
		}
	} 
	else if (change_type==6){
		// DOP Type 3: Single time-dependent DOPs obtained by adding fitness terms according to a set of templates	
		// Eq. 49			
		df_range=f_range*rho;
		for (k=0;k<ns;k++){			
			rand_perm(temp_v,perm_v,l);
			i=0;
			while (i<l){		
				if(i<os)					
					S_t30[perm_v[i]][k]=random_int(0,1);
				else 					
					S_t30[perm_v[i]][k]=-1;				// -1 indicates don´t care bits
				i++;		
			}
			
			df_t30[k]=random_dou()*(2*df_range)-df_range;			// observation: here, uniform distribution is used instead of normal distribution
		}
	}
	
	
	delete [] temp_v;
	delete [] perm_v;
}


