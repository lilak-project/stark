double calculate_solid_angle(double theta1, double theta2)
{
    double solid_angle = abs( 2*TMath::Pi() * (TMath::Cos(theta1*TMath::DegToRad())-TMath::Cos(theta2*TMath::DegToRad())) );
    return 0.8*solid_angle;
}

double calculate_yield(int days, double theta1, double theta2, double cross_section, double print_info=false)
{
    double solid_angle = calculate_solid_angle(theta1, theta2);
    double hours_a_day = 8;
    double total_sec = hours_a_day*60*60*days;
    double total_nbeam = 10E+5*total_sec;
    double td = 2.58E-8;
    double value = total_nbeam * td * solid_angle * cross_section;
    if (print_info)
        cout << "day=" << days 
             << ", yield ="  << setw(8) << value
             //<< ", tta ="    << setw(10) << Form("(%02d,%02d)",int(theta1),int(theta2))
             //<< ", n_beam =" << setw(10) << total_nbeam 
             //<< ", td ="     << setw(9) << td  
             //<< ", sa ="     << setw(9) << solid_angle 
             //<< ", cs ="     << setw(4) << cross_section 
             << ", stat_err = " << setw(4) << sqrt(value)/value*100
             << endl;
    return value;
}

void yield()
{
    //double d = 200; // mm
    //double d = 160; // mm
    double d = 160; // mm
    double r1 = 25.5; // mm
    double r2 = 81.95; // mm

    int tlab1 = atan2(r1,d) * TMath::RadToDeg();
    int tlab2 = atan2(r2,d) * TMath::RadToDeg();
    double tcom1 = 180 - 2 * tlab1; //atan2(r1,d) * TMath::RadToDeg();
    double tcom2 = 180 - 2 * tlab2; //atan2(r2,d) * TMath::RadToDeg();
    cout << "lab = " << tlab1 << " " << tlab2 << endl;
    cout << "com = " << tcom1 << " " << tcom2 << endl;

    //calculate_yield(1, 40,80,0.3);
    //calculate_yield(1, 68,72,0.3);
    //cout << endl;
    //calculate_yield(1, 100,120,7);
    //calculate_yield(1, 100,104,7);
    //cout << endl;
    //calculate_yield(1, tlab1,tlab2,3);
    //cout << endl;

    //calculate_yield(7, 68,72,0.3);
    //calculate_yield(7, 100,104,7);
    //calculate_yield(7, tlab1,tlab2,3);

    //cout << endl;
    //calculate_yield(1, 40,80,5);
    //calculate_yield(1, 100,120,16);
    //calculate_yield(1, tlab1,tlab2,0.25);
    //calculate_yield(1, tlab1,tlab1+4,0.25);

    //calculate_yield(1, 68,  72, 0.26);
    //calculate_yield(1, 100,104, 8.3);
    //calculate_yield(1, tlab1,tlab1+4, 1.9);

    //cout << endl;
    //calculate_yield(1, 68,  72, 5);
    //calculate_yield(1, 100,104, 6);
    //calculate_yield(1, tlab1,tlab1+4, 0.23);

    //cout << endl;
    //calculate_yield(1, 68,  72,       6);
    //calculate_yield(1, 100,104,       8.3);
    //calculate_yield(1, tlab1,tlab1+4, 3.1);

    //cout << endl;
    //calculate_yield(1, 68,  72,       7.7);
    //calculate_yield(1, 100,104,       6);
    //calculate_yield(1, tlab1,tlab1+4, 1.2);

    //double cs = 5.5;
    double cs = 1;
    cout << cs << endl;
    cout << endl;
    calculate_yield(4, 68,  72,       cs);
    calculate_yield(4, 100,104,       cs);
    calculate_yield(4, tlab1,tlab1+4, cs);

    cout << endl;
    calculate_yield(4, 68,  72,       cs);
    calculate_yield(4, 100,104,       cs);
    calculate_yield(4, tlab1,tlab1+4, cs);
}
