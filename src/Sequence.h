// Sequence class definition

class Sequence
{
  private:
         std::string id;
         std::string aa;
         int len;
         int aa_comp[23];
         double ext;
         double mol_wt;
         double mol_ext;
         double iso_point;
         double aliphatic;
         int atom_comp[5];
         double charge;
         int surfarea;
         double volume;
         double hydro;
         int num_atoms;
         int nneg;
         int npos;
         bool met;
         
   public:
         //constructors and defaults
         Sequence(){
                    for(unsigned int i=0; i < 23; i++)
                     aa_comp[i]=0;

                     hydro=volume=charge=iso_point=mol_wt=mol_ext=aliphatic=0.0;

                     num_atoms=surfarea=0;

                     for(unsigned int i=0; i<5; i++)
                      atom_comp[i] = 0;                    

                     met = 0;
                    };

         Sequence(const std::string& i, const std::string& a) : id(i), aa(a){  };

         //methods
         void set_properties();

         void set_composition();

         void calc_charge(double*, double&, double&);

         void calc_protons(double*, const double&);

         void calc_IEP();

         void set_id(const std::string& i){id = i;};

         void set_aa(const std::string& a){aa = a;};

         void print();

         

};
