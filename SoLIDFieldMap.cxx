//SoLIDTracking
#include "SoLIDFieldMap.h"

SoLIDFieldMap * SoLIDFieldMap::fInstance = NULL;

//__________________________________________________________________
SoLIDFieldMap::SoLIDFieldMap()
{
  for (int i=0; i<ZSIZE; i++){
    for (int j=0; j<RSIZE; j++) { 
      Bz[i][j] = 0;
      Br[i][j] = 0;
    }
  }
  fField.SetXYZ(0.,0.,0.);
  
  LoadFieldMap();
}
//__________________________________________________________________
void SoLIDFieldMap::LoadFieldMap()
{
  cout<<"loading SoLID field map"<<endl;
  ifstream infile;
  infile.open("solenoid_CLEOv8.dat");
  if (!infile.is_open()){
    cout<<"cannot open field map file"<<endl;
    exit(0);
  }
  double input[4];
  while(1){
    infile>>input[0]>>input[1]>>input[2]>>input[3];
    if (infile.eof()) break;
    
    int z_pos = fabs( (int)(input[1] + ZSHIFT) );
    int r_pos = fabs( (int)(input[0]) );
    Bz[z_pos][r_pos] = input[3]/1000.;//convert from gauss to kgauss
    Br[z_pos][r_pos] = input[2]/1000.;//convert from gauss to kgauss
  }
  
  infile.close();
  
}
//___________________________________________________________________
TVector3 & SoLIDFieldMap::GetBField(double x, double y, double z)
{
  //here use cm, other place use m
  x = 100*x;
  y = 100*y;
  z = 100*z + ZSHIFT;
  double r = sqrt(x*x + y*y);
  if (r >= RSIZE || z <= 0 || z >= ZSIZE){
    fField.Clear();
    return fField;
  }else{
  
    int z_max, z_min, r_max, r_min;
    r_min = (int)r;
    r_max = r_min + RSTEP;
    z_min = (int)z;
    z_max = z_min + ZSTEP;

    double f21_Bz =  Bz[z_max][r_min];
    double f21_Br =  Br[z_max][r_min];
    double f22_Bz =  Bz[z_max][r_max];
    double f22_Br =  Br[z_max][r_max];

    double f11_Bz =  Bz[z_min][r_min];
    double f11_Br =  Br[z_min][r_min];
    double f12_Bz =  Bz[z_min][r_max];
    double f12_Br =  Br[z_min][r_max];
    
    //linear interpolation
    double Bzi = (1./((r_max - r_min)*(z_max - z_min)))*(f11_Bz*( z_max - z )*( r_max - r ) +
                                                         f21_Bz*( z - z_min )*( r_max - r ) +
                                                         f12_Bz*( z_max - z )*( r - r_min ) +
                                                         f22_Bz*( z - z_min )*( r - r_min ) );

    double Bri = (1./((r_max - r_min)*(z_max - z_min)))*(f11_Br*( z_max - z )*( r_max - r ) +
                                                         f21_Br*( z - z_min )*( r_max - r ) +
                                                         f12_Br*( z_max - z )*( r - r_min ) +
                                                         f22_Br*( z - z_min )*( r - r_min ) );

   fField.SetXYZ(Bri*x/r, Bri*y/r, Bzi);
   return fField;
  }

}














