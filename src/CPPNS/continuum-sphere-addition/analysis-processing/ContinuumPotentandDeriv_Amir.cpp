double derivative (Frame &qqq, double &r, int &i){

   ofstream print ("derivative.log", ios::app);

   double R,r2, dx, dy, dz, sigma, epsilon, sigma2,rCutOff2, dudr, uR, density;
   
      if (i%10 == 0) {
    density =  0.298428946;  

    dudr = 0.0;
  
        //Calculating gas-solid interactions    
        for (int i = 0; i < qqq.nsolid ; i++)
            {
 
             for (int j = qqq.nsolid; j < qqq.natoms; j++)
               { 
                    dx = qqq.atom[i].x - qqq.atom[j].x;
                 dy = qqq.atom[i].y - qqq.atom[j].y;
                 dz = qqq.atom[i].z - qqq.atom[j].z;
                 pbc(dx, dy, dz, qqq.boxx, qqq.boxy, qqq.boxz);
                     r2 = pow (dx , 2.0) + pow (dy, 2.0) + pow (dz , 2.0);
                 R = sqrt (r2);
                      sigma = (qqq.atom[i].sigma + qqq.atom[j].sigma) / 2.0;
                      sigma2 = sigma * sigma;
                      epsilon = sqrt (qqq.atom[i].epsilon * qqq.atom[j].epsilon);

                    dudr += (16.0*density*epsilon*M_PI*pow(r,2)*pow(sigma,6)*(-5.0*pow(pow(r,2) - pow(R,2),6)*(pow(r,2) + pow(R,2)) +
       (5.0*pow(r,4) + 10.0*pow(r,2)*pow(R,2) + pow(R,4))*(pow(r,4) + 10.0*pow(r,2)*pow(R,2) + 5.0*pow(R,4))*pow(sigma,6)))/(5.0*pow(pow(r,2) - pow(R,2),10));
                   
              }
         }
        
         print << dudr << endl;
         }
}



double PotentialEnergy (Frame &qqq, double &a, double &density){
   int n;
   n = qqq.natoms;
   double R,r2, dx, dy, dz, sigma, epsilon, sigma2,rCutOff2, u, uR;
   
   
 
    uR = 4.0  * (pow ((1.0/qqq.rCutOff),12.0) - pow ((1.0/qqq.rCutOff), 6.0) );
    u = 0.0;
  
        //Calculating gas-solid interactions    
        for (int i = 0; i < qqq.nsolid ; i++)
            {
 
             for (int j = qqq.nsolid; j < qqq.natoms; j++)
               { 
                    dx = qqq.atom[i].x - qqq.atom[j].x;
                 dy = qqq.atom[i].y - qqq.atom[j].y;
                 dz = qqq.atom[i].z - qqq.atom[j].z;
                 pbc(dx, dy, dz, qqq.boxx, qqq.boxy, qqq.boxz);
                     r2 = pow (dx , 2.0) + pow (dy, 2.0) + pow (dz , 2.0);
                 R = sqrt (r2);
                      sigma = (qqq.atom[i].sigma + qqq.atom[j].sigma) / 2.0;
                      sigma2 = sigma * sigma;
                      epsilon = sqrt (qqq.atom[i].epsilon * qqq.atom[j].epsilon);

                    u += (16.0*pow(a,3)*density*epsilon*M_PI*pow(sigma,6)*(15.0*pow(pow(a,2) - pow(R,2),6) - (5.0*pow(a,6) + 45.0*pow(a,4)*pow(R,2) + 63.0*pow(a,2)*pow(R,4) + 15.0*pow(R,6))*pow(sigma,6)))/(45.0*pow(pow(a,2) - pow(R,2),9));
                   
              }
         }       
