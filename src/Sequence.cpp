// Sequence methods

#include <stdio.h>
#include <string>
#include <math.h>
#include "Sequence.h"


int atomC(const int& ch)
{
   static int carbons [] =
  {
      0, 3, 0, 3, 4, 5, 9, 2, 6, 6, 0, 6,
      6, 5, 4, 0, 5, 5, 6, 3, 4, 0, 5, 11,
      0, 9, 0
  };

 //   std::cout<<carbons[ch&31]<<std::endl;
    return (isalpha(ch) ? carbons[ch & 31] : 0 );
}

int atomN(const int& ch)
{
     static int nitrogens [] =
  {
      0, 1, 1, 1, 1, 1, 1, 1, 3, 1, 0, 2,
      1, 1, 2, 0, 1, 2, 4, 1, 1, 0, 1, 2,
      0, 1, 0
  };

  return (isalpha(ch) ? nitrogens[ch & 31] : 0 );
}

int atomO(const int& ch)
{
     static int oxygens [] =
  {
      0, 1, 0, 1, 3, 3, 1, 1, 1, 1, 0, 1, 1,
      1, 2, 0, 1, 2, 1, 2, 2, 0, 1, 1, 0, 2, 0
  };

  return (isalpha(ch) ? oxygens[ch & 31] : 0 );
}

int atomS(const int& ch)
{
     static int sulphurs [] =
  {
      0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0
  };
  return (isalpha(ch) ? sulphurs[ch & 31] : 0 );
}

int atomH(const int& ch)
{
     static int hydrogens [] =
  {
      0, 5, 0, 5, 5, 7, 9, 3, 7, 11, 0,
      12, 11, 9, 6, 0, 7, 8, 12, 5, 7, 0, 9,
      10, 0, 9, 0
  };

   return (isalpha(ch) ? hydrogens[ch & 31] : 0 );
}


int molExt(const int& ch)
{
   static int molEx [] =
  {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5690,
      0, 1280, 0
  };

    return (isalpha(ch) ? molEx[ch & 31] : 0 );

}

double pI(const int& ch)
{
   static double pIcvs [] =
  {
      0.0, 6.01, 0.0, 5.07, 2.77, 3.22, 5.48, 5.97, 7.59, 6.02, 0.0, 9.74,
      5.98, 5.74, 5.41, 0.0,  6.3, 5.65, 10.76, 5.68, 5.6, 0.0, 5.96, 5.89,
      0.0, 5.66, 0.0
  };

    return (isalpha(ch) ? pIcvs[ch & 31] : 0.0);
}

double Kd(const int& ch)
{
   static double dissoc [] =
  {
      0.0, 0.0, 0.0, -8.5, -3.9, -4.1, 0.0, 0.0, -6.5, 0.0, 0.0, -10.8,
      0.0, 0.0, 0.0, 0.0,  0.0, 0.0, -12.5, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, -10.1, 0.0
  };

  return (isalpha(ch) &&  dissoc[ch & 31] < 0.0 ? pow(10.0,dissoc[ch & 31]) : 0.0);
}
double pK(const int& ch)
{
  static double chargecvs [] =
  {
      0.0, 0.0, 0.0, 8.5, 3.9, 4.1, 0.0, 0.0, 6.5, 0.0, 0.0, 10.8,
      0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 12.5, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 10.1, 0.0
  };

    return (isalpha(ch) ? chargecvs[ch & 31] : 0.0);

}

int area(const int& ch)
{
    static int areacvs[] =
    {
      999, 115, 999, 135, 160, 180, 210,  75, 195, 175, 999, 200, 170, 185, 150,
      999, 145, 190, 225, 115, 140, 999, 140, 255, 999, 999, 999
    };
     
    return (isalpha(ch) ? areacvs[ch & 31] : 0);
}

double mwt(const int& ch)
{
    static double mwtcvs[] =
    {
      0, 89.09, 132.61, 121.12, 133.10, 147.13, 165.19, 75.07, 155.16, 131.17,
      0, 146.19, 131.17, 149.21, 132.12, 0, 115.13, 146.15, 174.2, 105.09,
      119.12, 0, 117.15, 204.23, 0, 181.19, 146.64
    };



    return (isalpha(ch) ? mwtcvs[ch & 31] : 0);
}

  



double kd_hydro(const int& ch)
{
  static double kd[] =
  {
      999.99, 1.8, 999.99, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, 0.00, -3.9, 3.8, 1.9, -3.5,
      0.00, -1.6, -3.5, -4.5, -0.8, -0.7, 0.00, 4.2, -0.9, 0.00, -1.3, 0.00
  };

    return (isalpha(ch) ? kd[ch & 31] : 0.00);
}

double vol(const int& ch)
{
    static double aavol[] =
    {
      0, 88.6, 999.99, 138.4, 114.1, 143.8, 189.9, 60.1, 153.2, 166.7, 999.99, 168.6, 166.7, 162.9, 111.1,
      0, 112.7, 108.5, 173.4, 89.0, 116.1, 0.0, 140.0, 227.8, 999.99, 193.6, 999.99
    };

    return (isalpha(ch) ? aavol[ch & 31] : 0.00 );
}

double pHtoConc(const double& pH)
{
  return pow(10.0,-pH);
}

int aanum(const int& ch)
{
    static int aacvs[] =
    {
      999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
      22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 22 );
}

int atoms(const int& ch)
{
    static int atomcvs[] =
    {
      999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
      22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };

    return (isalpha(ch) ? atomcvs[ch & 31] : 0);
}

void Sequence::calc_charge( double* p, double& c, double& h)
{
     double pamino   =  h/(h+pow(10,-8.6));
     double pcarboxy  =  h/(h+pow(10,-3.6));

     
     
     c = (p[1]+p[11]+p[8] + pamino) -
             ((double)aa_comp[18]     -  p[18]  +
              (double)aa_comp[4]      -  p[4]   +
              (double)aa_comp[3]      -  p[3]   +
              (double)aa_comp[6]      -  p[6]   +
              1.0 - pcarboxy );
    
}

void Sequence::calc_protons(double* p, const double& h)
{
  //get proton scale for amino acids at pH h
  std::string amino_acids = "ARNDCQEGHILKMFPSTWYVBZX";
  
  for(unsigned int i=0; i < amino_acids.length(); i++)
  {
    double kd = Kd(amino_acids[i]);
    //std::cout<<"Kd:  "<<kd<<std::endl;
    if(kd != 0.0)
    {
     p[i] =  (h/(h+Kd(amino_acids[i]))) * (double)aa_comp[aanum(amino_acids[i])];
     
    }
    else {
           p[i] = 0.0;
         }
  }
}

void Sequence::calc_IEP()
{
  double tpH = 1.0;
  double bpH = 14.0;
  
  double H   = pHtoConc(tpH);
  double* protons = new double[26];
  double top,btm,ch;

  calc_protons(protons,H);
  calc_charge(protons,top,H);

  H=pHtoConc(6.5);
  calc_protons(protons,H);
  calc_charge(protons,charge,H);

  H=pHtoConc(bpH);
  calc_protons(protons,H);
  calc_charge(protons,btm,H);
  

  double mid = 0.0;
  iso_point = tpH;
  
  while(bpH-tpH > 0.0001)
  {
     mid = ((bpH-tpH)/2.0) + tpH;

     H   = pHtoConc(mid);
     calc_protons(protons,H);
     calc_charge(protons,ch,H);

     

     if( ch > 0.0 )
     {
       tpH = mid;
       continue;
     }

     if( ch < 0.0 )
     {
       bpH = mid;
       continue;
     }
  }

  iso_point = mid;
   
}

void Sequence::set_composition()
{
   for(unsigned int i=0; i< aa.length(); i++)
   {
      
      aa_comp[aanum(aa[i])]++;
      
      num_atoms+= atoms(aa[i]);
      mol_ext  += molExt(aa[i]);
      mol_wt   += mwt(aa[i]);
      volume   += vol(aa[i]);
      surfarea += area(aa[i]);
      hydro    += kd_hydro(aa[i]);
      iso_point+= pI(aa[i]);

      atom_comp[0] += atomC(aa[i]);
      atom_comp[1] += atomH(aa[i]);
      atom_comp[2] += atomO(aa[i]);
      atom_comp[3] += atomN(aa[i]);
      atom_comp[4] += atomS(aa[i]);

     

   
   }

   //add on N and C terminal atoms H and OH
   atom_comp[1]+=2;
   atom_comp[2]++;
    
   //water correction to mol_wt
   mol_wt -= (18.015*(aa.length()-1));
   calc_IEP();

   num_atoms = atom_comp[0] + atom_comp[1] +
               atom_comp[2] + atom_comp[3] +
               atom_comp[4] ;
   double len= aa.length();   
   aliphatic =  aa_comp[0]          +
               ( 2.9 * aa_comp[19]) +
               ( 3.9 * (aa_comp[10] + aa_comp[9] ) );
   aliphatic *= 100/len;

   npos = aa_comp[1] + aa_comp[11];
   nneg = aa_comp[3] + aa_comp[6];
 
}


void Sequence::set_properties()
{
  if( aa[0] == 'M')
    met = true;
    
  set_composition();

}       

void Sequence::print()
{
  std::cout<<id<<"\t"
           <<met<<"\t"
           <<aa.length()<<"\t"
           <<mol_wt<<"\t"
           <<volume<<"\t"
           <<surfarea<<"\t"
           <<hydro<<"\t"
           <<hydro/(double)aa.length()<<"\t"
           <<charge<<"\t"
           <<mol_ext<<"\t"
           <<iso_point<<"\t"
           <<aliphatic<<"\t"
           <<npos<<"\t"
           <<nneg<<"\t"
           <<npos/(double)aa.length()<<"\t"
           <<nneg/(double)aa.length()<<"\t";

  for(unsigned int i=0; i< 23; i++)
  {
    std::cout<<aa_comp[i]<<"\t"<<aa_comp[i]/double(aa.length())<<"\t";
  }

  for(unsigned int i=0; i<5; i++)
  {
    std::cout<<atom_comp[i]<<"\t"<<atom_comp[i]/double(num_atoms)<<"\t";
  }
  std::cout<<num_atoms<<std::endl;
  
}
