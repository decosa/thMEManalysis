//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.4.8, 2012-07-24
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "TopQuarkAnalysis/SingleTop/interface/HelAmps_sm_no_b_mass.h"

namespace MG5_sm_no_b_mass 
{


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
}

void ixxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fi[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 
  fi[4] = complex<double> (p[0] * nsf, p[3] * nsf); 
  fi[5] = complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.0)
  {
    pp = min(p[0], pow((pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2)), 0.5)); 
    if (pp == 0.0)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      fi[0] = ip * sqm[ip]; 
      fi[1] = im * nsf * sqm[ip]; 
      fi[2] = ip * nsf * sqm[im]; 
      fi[3] = im * sqm[im]; 
    }
    else
    {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.0); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fi[0] = sfomega[0] * chi[im]; 
      fi[1] = sfomega[0] * chi[ip]; 
      fi[2] = sfomega[1] * chi[im]; 
      fi[3] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.0), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = complex<double> (-nhel * pow(2.0 * p[0], 0.5), 0.0); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[0] = complex<double> (0.0, 0.0); 
      fi[1] = complex<double> (0.0, 0.0); 
      fi[2] = chi[0]; 
      fi[3] = chi[1]; 
    }
    else
    {
      fi[0] = chi[1]; 
      fi[1] = chi[0]; 
      fi[2] = complex<double> (0.0, 0.0); 
      fi[3] = complex<double> (0.0, 0.0); 
    }
  }
  return; 
}

void sxxxxx(double p[4], int nss, complex<double> sc[3])
{
  sc[0] = complex<double> (1.00, 0.00); 
  sc[1] = complex<double> (p[0] * nss, p[3] * nss); 
  sc[2] = complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}

void oxxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fo[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 
  fo[4] = complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[5] = complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
    if (pp == 0.000)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = -((1 - nh)/2) * nhel; 
      im = (1 + nh)/2 * nhel; 
      fo[0] = im * sqm[im]; 
      fo[1] = ip * nsf * sqm[im]; 
      fo[2] = im * nsf * sqm[ - ip]; 
      fo[3] = ip * sqm[ - ip]; 
    }
    else
    {
      pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.00); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], -p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fo[0] = sfomeg[1] * chi[im]; 
      fo[1] = sfomeg[1] * chi[ip]; 
      fo[2] = sfomeg[0] * chi[im]; 
      fo[3] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.00), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.00); 
    if(sqp0p3 == 0.000)
    {
      chi[1] = complex<double> (-nhel, 0.00) * pow(2.0 * p[0], 0.5); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fo[0] = chi[0]; 
      fo[1] = chi[1]; 
      fo[2] = complex<double> (0.00, 0.00); 
      fo[3] = complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[0] = complex<double> (0.00, 0.00); 
      fo[1] = complex<double> (0.00, 0.00); 
      fo[2] = chi[1]; 
      fo[3] = chi[0]; 
    }
  }
  return; 
}

void vxxxxx(double p[4], double vmass, int nhel, int nsv, complex<double> vc[6])
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 
  sqh = pow(0.5, 0.5); 
  hel = double(nhel); 
  nsvahl = nsv * abs(hel); 
  pt2 = pow(p[1], 2) + pow(p[2], 2); 
  pp = min(p[0], pow(pt2 + pow(p[3], 2), 0.5)); 
  pt = min(pp, pow(pt2, 0.5)); 
  vc[4] = complex<double> (p[0] * nsv, p[3] * nsv); 
  vc[5] = complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - abs(hel); 
    if(pp == 0.0)
    {
      vc[0] = complex<double> (0.0, 0.0); 
      vc[1] = complex<double> (-hel * sqh, 0.0); 
      vc[2] = complex<double> (0.0, nsvahl * sqh); 
      vc[3] = complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[0] = complex<double> (hel0 * pp/vmass, 0.0); 
      vc[3] = complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[1] = complex<double> (hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
            p[2]/pt * sqh);
        vc[2] = complex<double> (hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
            p[1]/pt * sqh);
      }
      else
      {
        vc[1] = complex<double> (-hel * sqh, 0.0); 
        vc[2] = complex<double> (0.0, nsvahl * Sgn(sqh, p[3])); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = pow(pow(p[1], 2) + pow(p[2], 2), 0.5); 
    vc[0] = complex<double> (0.0, 0.0); 
    vc[3] = complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[1] = complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[2] = complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[1] = complex<double> (-hel * sqh, 0.0); 
      vc[2] = complex<double> (0.0, nsv * Sgn(sqh, p[3])); 
    }
  }
  return; 
}

void FFS1_1(complex<double> F2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> denom; 
  double P1[4]; 
  F1[4] = F2[4] + S3[1]; 
  F1[5] = F2[5] + S3[2]; 
  P1[0] = F1[4].real(); 
  P1[1] = F1[5].real(); 
  P1[2] = F1[5].imag(); 
  P1[3] = F1[4].imag(); 
  denom = 1./(((M1 * (-M1 + complex<double> (0., 1.) * W1)) + ((pow(P1[0], 2))
      - (pow(P1[1], 2)) - (pow(P1[2], 2)) - (pow(P1[3], 2)))));
  F1[0] = COUP * denom * (S3[0] * ((F2[3] * (complex<double> (0., 1.) * P1[1] -
      P1[2])) + ((F2[2] * (complex<double> (0., 1.) * P1[0] + complex<double>
      (0., 1.) * P1[3])) + complex<double> (0., 1.) * (F2[0] * M1))));
  F1[1] = COUP * denom * (S3[0] * ((F2[2] * (complex<double> (0., 1.) * P1[1] +
      P1[2])) + ((F2[3] * (complex<double> (0., 1.) * P1[0] + complex<double>
      (0., -1.) * P1[3])) + complex<double> (0., 1.) * (F2[1] * M1))));
  F1[2] = COUP * denom * (S3[0] * ((F2[1] * (complex<double> (0., -1.) * P1[1]
      + P1[2])) + ((F2[0] * (complex<double> (0., 1.) * P1[0] + complex<double>
      (0., -1.) * P1[3])) + complex<double> (0., 1.) * (F2[2] * M1))));
  F1[3] = COUP * denom * (S3[0] * ((F2[0] * (complex<double> (0., -1.) * P1[1]
      - P1[2])) + ((F2[1] * (complex<double> (0., 1.) * P1[0] + complex<double>
      (0., 1.) * P1[3])) + complex<double> (0., 1.) * (F2[3] * M1))));
}


void FFV2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> denom; 
  complex<double> OM3; 
  double P3[4]; 
  V3[4] = -F1[4] + F2[4]; 
  V3[5] = -F1[5] + F2[5]; 
  P3[0] = -V3[4].real(); 
  P3[1] = -V3[5].real(); 
  P3[2] = -V3[5].imag(); 
  P3[3] = -V3[4].imag(); 
  OM3 = 0; 
  if (M3 != 0)
    OM3 = 1./pow(M3, 2); 
  denom = 1./(((M3 * (-M3 + complex<double> (0., 1.) * W3)) + ((pow(P3[0], 2))
      - (pow(P3[1], 2)) - (pow(P3[2], 2)) - (pow(P3[3], 2)))));
  V3[0] = COUP * denom * ((OM3 * ((F1[0] * ((F2[3] * (complex<double> (0., 1.)
      * P3[1] - P3[2])) + (F2[2] * (complex<double> (0., 1.) * P3[0] +
      complex<double> (0., 1.) * P3[3])))) + (F1[1] * ((F2[2] *
      (complex<double> (0., 1.) * P3[1] + P3[2])) + (F2[3] * (complex<double>
      (0., 1.) * P3[0] + complex<double> (0., -1.) * P3[3]))))) * P3[0]) +
      (complex<double> (0., -1.) * (F2[2] * F1[0]) + complex<double> (0., -1.)
      * (F2[3] * F1[1])));
  V3[1] = COUP * denom * ((OM3 * ((F1[0] * ((F2[3] * (complex<double> (0., 1.)
      * P3[1] - P3[2])) + (F2[2] * (complex<double> (0., 1.) * P3[0] +
      complex<double> (0., 1.) * P3[3])))) + (F1[1] * ((F2[2] *
      (complex<double> (0., 1.) * P3[1] + P3[2])) + (F2[3] * (complex<double>
      (0., 1.) * P3[0] + complex<double> (0., -1.) * P3[3]))))) * P3[1]) +
      (complex<double> (0., 1.) * (F2[3] * F1[0]) + complex<double> (0., 1.) *
      (F2[2] * F1[1])));
  V3[2] = COUP * denom * ((OM3 * ((F1[0] * ((F2[3] * (complex<double> (0., 1.)
      * P3[1] - P3[2])) + (F2[2] * (complex<double> (0., 1.) * P3[0] +
      complex<double> (0., 1.) * P3[3])))) + (F1[1] * ((F2[2] *
      (complex<double> (0., 1.) * P3[1] + P3[2])) + (F2[3] * (complex<double>
      (0., 1.) * P3[0] + complex<double> (0., -1.) * P3[3]))))) * P3[2]) +
      (-(F2[3] * F1[0]) + (F2[2] * F1[1])));
  V3[3] = COUP * denom * ((OM3 * ((F1[0] * ((F2[3] * (complex<double> (0., 1.)
      * P3[1] - P3[2])) + (F2[2] * (complex<double> (0., 1.) * P3[0] +
      complex<double> (0., 1.) * P3[3])))) + (F1[1] * ((F2[2] *
      (complex<double> (0., 1.) * P3[1] + P3[2])) + (F2[3] * (complex<double>
      (0., 1.) * P3[0] + complex<double> (0., -1.) * P3[3]))))) * P3[3]) +
      (complex<double> (0., 1.) * (F2[2] * F1[0]) + complex<double> (0., -1.) *
      (F2[3] * F1[1])));
}


void FFV2_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  vertex = COUP * ((F2[2] * ((F1[0] * (complex<double> (0., -1.) * V3[0] +
      complex<double> (0., -1.) * V3[3])) + (F1[1] * (complex<double> (0., -1.)
      * V3[1] - V3[2])))) + (F2[3] * ((F1[0] * (complex<double> (0., -1.) *
      V3[1] + V3[2])) + (F1[1] * (complex<double> (0., -1.) * V3[0] +
      complex<double> (0., 1.) * V3[3])))));
}


void FFS1_2(complex<double> F1[], complex<double> S3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> denom; 
  double P2[4]; 
  F2[4] = F1[4] - S3[1]; 
  F2[5] = F1[5] - S3[2]; 
  P2[0] = F2[4].real(); 
  P2[1] = F2[5].real(); 
  P2[2] = F2[5].imag(); 
  P2[3] = F2[4].imag(); 
  denom = 1./(((M2 * (-M2 + complex<double> (0., 1.) * W2)) + ((pow(P2[0], 2))
      - (pow(P2[1], 2)) - (pow(P2[2], 2)) - (pow(P2[3], 2)))));
  F2[0] = COUP * denom * (S3[0] * ((F1[3] * (complex<double> (0., -1.) * P2[1]
      - P2[2])) + ((F1[2] * (complex<double> (0., 1.) * P2[0] + complex<double>
      (0., -1.) * P2[3])) + complex<double> (0., 1.) * (F1[0] * M2))));
  F2[1] = COUP * denom * (S3[0] * ((F1[2] * (complex<double> (0., -1.) * P2[1]
      + P2[2])) + ((F1[3] * (complex<double> (0., 1.) * P2[0] + complex<double>
      (0., 1.) * P2[3])) + complex<double> (0., 1.) * (F1[1] * M2))));
  F2[2] = COUP * denom * (S3[0] * ((F1[1] * (complex<double> (0., 1.) * P2[1] +
      P2[2])) + ((F1[0] * (complex<double> (0., 1.) * P2[0] + complex<double>
      (0., 1.) * P2[3])) + complex<double> (0., 1.) * (F1[2] * M2))));
  F2[3] = COUP * denom * (S3[0] * ((F1[0] * (complex<double> (0., 1.) * P2[1] -
      P2[2])) + ((F1[1] * (complex<double> (0., 1.) * P2[0] + complex<double>
      (0., -1.) * P2[3])) + complex<double> (0., 1.) * (F1[3] * M2))));
}


void VVS1_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  vertex = COUP * (S3[0] * (complex<double> (0., -1.) * (V2[0] * V1[0]) +
      complex<double> (0., 1.) * (V2[1] * V1[1]) + complex<double> (0., 1.) *
      (V2[2] * V1[2]) + complex<double> (0., 1.) * (V2[3] * V1[3])));
}


}  // end namespace $(namespace)s_sm_no_b_mas

