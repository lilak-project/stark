void print_data()
{
     for (int energy : {4,5,8})
     {
        TString fn = Form("data_x/cross_section_scaled_e%d.txt",energy);
        ifstream file(fn);
        double x, y, ex, ey;
        cout << endl;
        int percentage_count = 0;
        double percentage_avg = 0;
        while (file >> x >> ex >> y >> ey) {
            auto percentage = ey/y*100;
            percentage_avg += percentage;
            cout << energy << setw(4) << x << setw(10) << Form("%.1f",percentage) << " %" << endl;
            percentage_count++;
        }
        percentage_avg = percentage_avg/percentage_count;
        cout << "(mean = " << Form("%.1f",percentage_avg) << " %)" << endl;
     }
}
