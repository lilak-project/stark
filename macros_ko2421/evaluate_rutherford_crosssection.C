#ifndef evaluate_rutherford_crosssection_C
#define evaluate_rutherford_crosssection_C

const double e_squared_MeVfm = 1.44; // in MeV·fm (Coulomb constant (k) * e^2)
const double fm2_to_mbarn = 10;
const double m_34Ar = 33.967867; // u
const double m_proton = 1.007276; // u
const int z_34Ar = 18;
const int z_proton = 1;

double EvalRutherford(double E_lab, double theta_com_deg)
{
    // rutherford cross section
    // = ( (e^2/(4*pi*epsilon)) * Z1*Z2 / (4*E_com) )**2 / sin^4(theta/2) [fm2]
    // = ( e_squared_MeVfm * zzo4e )**2 / sinp4 [fm2]
    // = ( e_squared_MeVfm * zzo4e )**2 / sinp4 * fm2_to_mbarn [mbarn/sr]
    double E_com = E_lab * (m_34Ar / (m_34Ar + m_proton));
    double theta = theta_com_deg * TMath::DegToRad();
    double zzo4e = (z_34Ar * z_proton) / (4. * E_com);
    double sinp4 = pow(sin(theta/2.), 4);
    double cross_section_fm2 = pow(zzo4e * e_squared_MeVfm, 2) / sinp4;
    double cross_section_mbarn = cross_section_fm2 * fm2_to_mbarn;
    return cross_section_mbarn; // in mbarn/sr
}

const double constant1 = 1.6021E-19; // 전하의 크기 (Coulombs)
const double constant2 = 8.854E-12; // 진공유전율
const double constant3 = 1.0E+31; // 1m2 = 1E+31 mbarn
const double constant4 = 1.602E-13; // MeV to J

double EvalRutherford_DT(double E_lab, double theta_com_deg)
{
    double factor = 1 * 18 * TMath::Power(constant1, 2) / constant2 / 4 / 4 / TMath::Pi() / (E_lab * constant4);
    double numerator = TMath::Power(factor, 2);
    double denominator = TMath::Power(TMath::Sin(TMath::DegToRad()*theta_com_deg / 2), 4);
    double rutherford = numerator / denominator * constant3;
    return rutherford;
}

#endif

//double EvalRutherford(double E_lab, double theta_com_deg)
//{
//    double theta_com_rad = theta_com_deg * TMath::Pi() / 180.0;
//    double sin_half_theta = sin(theta_com_rad / 2.0);
//    double E_com = E_lab * (m_34Ar / (m_34Ar + m_proton));
//    double factor = (z_34Ar * z_proton * eSquared) / (4.0 * E_com);
//    double cross_section_fm2 = pow(factor, 2) / pow(sin_half_theta, 4);
//    double cross_section_mbarn = cross_section_fm2 * fm2_to_mbarn;
//    return cross_section_mbarn; // in mbarn/sr
//}
