void draw_summary()
{
    auto cvs = new TCanvas("cvs","",2000*1,1200*1);
    cvs -> Divide(3,2);

    auto top4 = new LKDrawingGroup("data_lilak/anaX_e4_cutX_ProtonCut.root",false);
    auto draw1 = top4 -> FindGroup("x_final") -> FindDrawing("drawX_99"); draw1 -> SetCanvas(cvs -> cd(1));  draw1 -> Draw();
    auto draw4 = top4 -> FindGroup("x_final") -> FindDrawing("drawR_99"); draw4 -> SetCanvas(cvs -> cd(4));  draw4 -> Draw();
    auto top5 = new LKDrawingGroup("data_lilak/anaX_e5_cutX_ProtonCut.root",false);
    auto draw2 = top5 -> FindGroup("x_final") -> FindDrawing("drawX_99"); draw2 -> SetCanvas(cvs -> cd(2));  draw2 -> Draw();
    auto draw5 = top5 -> FindGroup("x_final") -> FindDrawing("drawR_99"); draw5 -> SetCanvas(cvs -> cd(5));  draw5 -> Draw();
    auto top8 = new LKDrawingGroup("data_lilak/anaX_e8_cutX_ProtonCut.root",false);
    auto draw3 = top8 -> FindGroup("x_final") -> FindDrawing("drawX_99"); draw3 -> SetCanvas(cvs -> cd(3));  draw3 -> Draw();
    auto draw6 = top8 -> FindGroup("x_final") -> FindDrawing("drawR_99"); draw6 -> SetCanvas(cvs -> cd(6));  draw6 -> Draw();

    cvs -> SaveAs("data_x/summary.png");
}
