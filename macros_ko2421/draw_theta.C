void draw_theta()
{
    auto file = new TFile("data_viewer/et_All.dummy.root");
    auto group = file -> Get("et_All");
    group -> Draw();
}
