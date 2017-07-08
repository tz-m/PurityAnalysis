#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TBox.h"
#include "TGaxis.h"
#include "TFrame.h"

void Put(std::map<Int_t,std::map<Int_t,Int_t> > & map, Int_t tpc, Int_t wire)
{
  if (map.find(tpc) == map.end())
    {
      map[tpc][wire] = 1;
    }
  else
    {
      if (map[tpc].find(wire) == map[tpc].end())
        {
          map[tpc][wire] = 1;
        }
      else
        {
          map[tpc][wire] += 1;
        }
    }
}

Int_t HasDuplicate(std::vector<Int_t> hitvec)
{
  std::vector<Int_t> hitveccopy(hitvec.begin(),hitvec.end());
  std::sort(hitveccopy.begin(),hitveccopy.end());
  Int_t numdup = 0;
  for (size_t i = 0; i < hitveccopy.size(); ++i)
    {
      if (i>0 && hitveccopy[i] == hitveccopy[i-1]) {
          std::cout << "    wire " << hitveccopy[i] << std::endl;
          ++numdup;
        }
    }
  return numdup;
}

void RecoEfficiency()
{
  //TFile * file = TFile::Open("/media/mthiesse/Dell Portable Hard Drive/PurityData/robustreco_data_hist.root","READ");
  TFile * file = TFile::Open("/media/mthiesse/Dell Portable Hard Drive/PurityData/robust_oldgain_3000us_mcscale1.2_hist2.root","READ");

  TTreeReader reader("robusthit/RobustHitFinder",file);
  TTreeReaderValue<std::vector<Bool_t> > ismctruth(reader,"ismctruth");
  TTreeReaderArray<Float_t> hitt(reader,"hitt");
  TTreeReaderArray<Float_t> integral(reader,"integral");
  TTreeReaderArray<Float_t> segmentlength(reader,"segmentlength");
  TTreeReaderValue<UInt_t> c1(reader,"c1");
  TTreeReaderValue<UInt_t> c2(reader,"c2");
  TTreeReaderValue<UInt_t> trignum(reader,"trignum");
  TTreeReaderValue<Int_t> run(reader,"run");
  TTreeReaderValue<Int_t> event(reader,"event");
  TTreeReaderValue<Bool_t> fitsuccess(reader,"fitsuccess");
  TTreeReaderValue<Float_t> fitchi2(reader,"fitchi2");
  TTreeReaderValue<Float_t> fitconstant(reader,"fitconstant");
  TTreeReaderValue<Float_t> fitndf(reader,"fitndf");
  TTreeReaderValue<Float_t> fitquadratic(reader,"fitquadratic");
  TTreeReaderValue<Float_t> fitsumsqrresidual(reader,"fitsumsqrresidual");
  TTreeReaderArray<Int_t> allchannels(reader,"allchannels");
  TTreeReaderArray<Int_t> channel(reader,"channel");
  TTreeReaderArray<Int_t> tpc(reader,"tpc");
  TTreeReaderArray<Int_t> wire(reader,"wire");
  TTreeReaderArray<Int_t> width(reader,"width");
  TTreeReaderArray<Float_t> rms(reader,"rms");
  TTreeReaderArray<Float_t> baseline(reader,"baseline");
  TTreeReaderValue<std::vector<Bool_t> > assumedhit(reader,"assumedhit");
  TTreeReaderValue<std::vector<Bool_t> > fitrealhit(reader,"fitrealhit");
  TTreeReaderValue<Int_t> nwiresTPC0(reader,"nwiresTPC0");
  TTreeReaderValue<Int_t> nwiresTPC1(reader,"nwiresTPC1");
  TTreeReaderValue<Int_t> nwiresTPC2(reader,"nwiresTPC2");
  TTreeReaderValue<Int_t> nwiresTPC3(reader,"nwiresTPC3");
  TTreeReaderValue<Int_t> nwiresTPC4(reader,"nwiresTPC4");
  TTreeReaderValue<Int_t> nwiresTPC5(reader,"nwiresTPC5");
  TTreeReaderValue<Int_t> nwiresTPC6(reader,"nwiresTPC6");
  TTreeReaderValue<Int_t> nwiresTPC7(reader,"nwiresTPC7");

  TH1::SetDefaultSumw2(kTRUE);

  TH1I * tpc1 = new TH1I("tpc1","TPC 1",112,-0.5,111.5);
  TH1I * tpc3 = new TH1I("tpc3","TPC 3",112,-0.5,111.5);
  TH1I * tpc5 = new TH1I("tpc5","TPC 5",112,-0.5,111.5);
  TH1I * tpc7 = new TH1I("tpc7","TPC 7",112,-0.5,111.5);

  TH1F * eff1 = new TH1F("eff1","eff 1",25,0,1);
  TH1F * eff3 = new TH1F("eff3","eff 3",25,0,1);
  TH1F * eff5 = new TH1F("eff5","eff 5",25,0,1);
  TH1F * eff7 = new TH1F("eff7","eff 7",25,0,1);

  std::map<Int_t,std::map<Int_t,Int_t> > tpc_wirehitsNnumerator;
  std::map<Int_t,std::map<Int_t,Int_t> > tpc_wirehitsNdenominator;

  TGraphErrors * eff1wires = new TGraphErrors();
  TGraphErrors * eff3wires = new TGraphErrors();
  TGraphErrors * eff5wires = new TGraphErrors();
  TGraphErrors * eff7wires = new TGraphErrors();

  std::vector<Int_t> badchannels = {2, 4, 9, 31, 36, 39, 51, 65, 80, 84, 88, 92, 92, 100, 101, 102, 103, 104, 105, 106, 111, 110, 117, 120, 145, 146, 147, 148, 183, 204, 213, 215, 216, 217, 218, 219, 238, 236, 240, 249, 252, 253, 254, 255, 256, 257, 258, 259, 260, 280, 287, 319, 346, 347, 348, 357, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 432, 433, 434, 441, 444, 446, 451, 458, 460, 469, 493, 496, 500, 513, 517, 520, 524, 530, 533, 538, 539, 547, 548, 560, 561, 564, 565, 566, 567, 568, 569, 560, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 611, 619, 620, 621, 622, 623, 626, 632, 649, 650, 651, 652, 653, 654, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 669, 675, 692, 693, 694, 728, 729, 730, 731, 732, 736, 738, 744, 760, 771, 776, 788, 792, 804, 809, 819, 821, 822, 826, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 826, 826, 826, 826, 826, 826, 826, 884, 885, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 940, 943, 948, 952, 965, 978, 980, 992, 993, 1016, 1021, 1023, 1025, 1026, 1047, 1051, 1055, 1058, 1062, 1063, 1077, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1097, 1102, 1105, 1122, 1126, 1134, 1141, 1145, 1149, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176, 1189, 1194, 1204, 1205, 1206, 1207, 1208, 1209, 1210, 1211, 1216, 1217, 1226, 1240, 1241, 1242, 1243, 1244, 1245, 1246, 1247, 1254, 1255, 1259, 1272, 1273, 1279, 1295, 1298, 1299, 1310, 1312, 1313, 1314, 1315, 1316, 1317, 1318, 1319, 1320, 1321, 1322, 1323, 1324, 1325, 1340, 1358, 1378, 1382, 1388, 1397, 1407, 1413, 1421, 1430, 1439, 1439, 1440, 1442, 1447, 1450, 1458, 1461, 1468, 1472, 1472, 1478, 1486, 1494, 1495, 1496, 1497, 1498, 1499, 1500, 1501, 1502, 1503, 1504, 1505, 1506, 1507, 1510, 1511, 1520, 1528, 1546, 1547, 1554, 1555, 1556, 1557, 1558, 1559, 1560, 1561, 1562, 1563, 1564, 1565, 1566, 1567, 1568, 1569, 1570, 1571, 1574, 1575, 1588, 1589, 1599, 1606, 1607, 1617, 1621, 1625, 1637, 1638, 1639, 1640, 1641, 1642, 1643, 1646, 1647, 1660, 1661, 1671, 1672, 1673, 1674, 1675, 1676, 1677, 1678, 1679, 1694, 1716, 1717, 1718, 1719, 1720, 1721, 1722, 1723, 1724, 1725, 1726, 1727, 1728, 1729, 1730, 1731, 1732, 1733, 1740, 1741, 1742, 1752, 1753, 1754, 1755, 1756, 1757, 1758, 1759, 1788, 1792, 1824, 1825, 1826, 1827, 1828, 1829, 1830, 1831, 1832, 1833, 1834, 1835, 1836, 1837, 1843, 1866, 1870, 1880, 1881, 1882, 1883, 1884, 1885, 1886, 1887, 1888, 1889, 1890, 1891, 1892, 1893, 1894, 1895, 1896, 1897, 1898, 1899, 1900, 1901, 1902, 1903, 1904, 1905, 1906, 1907, 1909, 1936, 1946, 1948, 1950, 1951, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1974, 1992, 1995, 2004, 2029, 2034, 2035, 2036, 2037, 2038, 2039, 2040, 2041, 2042, 2043, 2044, 2045, 2046, 2047};

  std::vector<Int_t> edgechannels = {288, 399, 400, 511, 800, 911, 912, 1023, 1312, 1423, 1424, 1535, 1824, 1935, 1936, 2047};

  std::vector<Int_t> nwires1;
  std::vector<Int_t> nwires3;
  std::vector<Int_t> nwires5;
  std::vector<Int_t> nwires7;

  std::vector<TBox*> box1;
  std::vector<TBox*> box3;
  std::vector<TBox*> box5;
  std::vector<TBox*> box7;
  std::vector<TBox*> nbox1;
  std::vector<TBox*> nbox3;
  std::vector<TBox*> nbox5;
  std::vector<TBox*> nbox7;

  Int_t startch[8] = {0,400,0,912,0,1424,0,1936};

  for (int itpc = 1; itpc<=7; itpc+=2)
    {
      for (int iw = 0; iw < 112; ++iw)
        {
          if ((std::find(badchannels.begin(),badchannels.end(),iw+startch[itpc]) != badchannels.end()) ||
              (std::find(edgechannels.begin(),edgechannels.end(),iw+startch[itpc]) != edgechannels.end()))
            {
              TBox * b = new TBox((Double_t)iw-0.5,0,(Double_t)iw+0.5,1.0);
              b->SetFillColor(kRed);
              b->SetFillColorAlpha(kRed,0.4);
              TBox * nb = new TBox((Double_t)iw-0.5,0,(Double_t)iw+0.5,20000);
              nb->SetFillColor(kRed);
              nb->SetFillColorAlpha(kRed,0.4);
              switch (itpc)
                {
                case 1: box1.push_back(b); nbox1.push_back(nb); break;
                case 3: box3.push_back(b); nbox3.push_back(nb); break;
                case 5: box5.push_back(b); nbox5.push_back(nb); break;
                case 7: box7.push_back(b); nbox7.push_back(nb); break;
                }
            }
        }
    }

  while (reader.Next())
    {
      if (*run > 17512) continue;
      if (*trignum != 111) continue;
      if ( (*c1 == 28 && (*c2 == 6 || *c2 == 7 || *c2 == 8 || *c2 == 9)) ||
           (*c1 == 29 && (*c2 == 6 || *c2 == 7 || *c2 == 8)) ||
           (*c1 == 30 && (*c2 == 6 || *c2 == 7)) ||
           (*c1 == 31 && (*c2 == 6)) ||
           (*c1 == 36 && (*c2 == 15)) ||
           (*c1 == 37 && (*c2 == 15 || *c2 == 14))) continue;
      if ( (*fitsumsqrresidual / *fitndf) > 3 ) continue;
      if (!(*fitsuccess)) continue;

      Int_t i = 0;
      Int_t tpc1minwire=999, tpc1maxwire=-999, tpc3minwire=999, tpc3maxwire=-999, tpc5minwire=999, tpc5maxwire=-999, tpc7minwire=999, tpc7maxwire=-999;
      std::vector<Int_t> tpc1hits;
      std::vector<Int_t> tpc3hits;
      std::vector<Int_t> tpc5hits;
      std::vector<Int_t> tpc7hits;
      std::vector<Int_t> currentchannels;
      if (std::find(nwires1.begin(),nwires1.end(),*nwiresTPC1) == nwires1.end()) nwires1.push_back(*nwiresTPC1);
      if (std::find(nwires3.begin(),nwires3.end(),*nwiresTPC3) == nwires3.end()) nwires3.push_back(*nwiresTPC3);
      if (std::find(nwires5.begin(),nwires5.end(),*nwiresTPC5) == nwires5.end()) nwires5.push_back(*nwiresTPC5);
      if (std::find(nwires7.begin(),nwires7.end(),*nwiresTPC7) == nwires7.end()) nwires7.push_back(*nwiresTPC7);
      for (auto cit = allchannels.begin(); cit != allchannels.end(); ++cit)
        {
          currentchannels.push_back(*cit);
        }
      for (auto wit = wire.begin(); wit != wire.end(); ++wit)
        {
          if ((*ismctruth)[i]) continue;
          int wirenum = *wit;
          if ((*fitrealhit)[i] || (*assumedhit)[i])
            {
              if (tpc[i] == 1) {
                  if (wirenum != 0 && wirenum != 111) tpc1->Fill(wirenum);
                  if (wirenum < tpc1minwire) tpc1minwire = wirenum;
                  if (wirenum > tpc1maxwire) tpc1maxwire = wirenum;
                  tpc1hits.push_back(wirenum);
                }
              else if (tpc[i] == 3) {
                  if (wirenum != 0 && wirenum != 111) tpc3->Fill(wirenum);
                  if (wirenum < tpc3minwire) tpc3minwire = wirenum;
                  if (wirenum > tpc3maxwire) tpc3maxwire = wirenum;
                  tpc3hits.push_back(wirenum);
                }
              else if (tpc[i] == 5) {
                  if (wirenum != 0 && wirenum != 111) tpc5->Fill(wirenum);
                  if (wirenum < tpc5minwire) tpc5minwire = wirenum;
                  if (wirenum > tpc5maxwire) tpc5maxwire = wirenum;
                  tpc5hits.push_back(wirenum);
                }
              else if (tpc[i] == 7) {
                  if (wirenum != 0 && wirenum != 111) tpc7->Fill(wirenum);
                  if (wirenum < tpc7minwire) tpc7minwire = wirenum;
                  if (wirenum > tpc7maxwire) tpc7maxwire = wirenum;
                  tpc7hits.push_back(wirenum);
                }
            }
          ++i;
        }
      std::string hits3string = "TPC3 >> ";
      std::string hits5string = "TPC5 >> ";
      Int_t count_hits1 = 0;
      Int_t count_wires1 = 0;
      //if (HasDuplicate(tpc1hits)>0) std::cout << " in TPC1" << std::endl;
      //if (HasDuplicate(tpc3hits)>0) std::cout << " in TPC3" << std::endl;
      //if (HasDuplicate(tpc5hits)>0) std::cout << " in TPC5" << std::endl;
      //if (HasDuplicate(tpc7hits)>0) std::cout << " in TPC7" << std::endl;

      if (*nwiresTPC1 != 0)
        {
          std::vector<Int_t> allwire_tpc1;
          std::vector<Int_t> hits_tpc1;
          for (int iw = 0; iw < 112; ++iw)
            {
              if ((std::find(badchannels.begin(),badchannels.end(),iw+400) != badchannels.end()) ||
                  (std::find(edgechannels.begin(),edgechannels.end(),iw+400) != edgechannels.end())) continue;
              if (std::find(currentchannels.begin(),currentchannels.end(),iw+400) == currentchannels.end()) continue;
              //if ((std::find(tpc1hits.begin(),tpc1hits.end(),iw) != tpc1hits.end()) && (std::find(currentchannels.begin(),currentchannels.end(),iw+400) == currentchannels.end())) continue;
              if (std::find(currentchannels.begin(),currentchannels.end(),iw+400) != currentchannels.end()) allwire_tpc1.push_back(iw);
              if (std::find(tpc1hits.begin(),tpc1hits.end(),iw) != tpc1hits.end()) hits_tpc1.push_back(iw);
              if (std::find(tpc1hits.begin(),tpc1hits.end(),iw) != tpc1hits.end()) {
                  ++count_hits1;
                  Put(tpc_wirehitsNnumerator,1,iw);
                }
              ++count_wires1;
              Put(tpc_wirehitsNdenominator,1,iw);
            }
          for (auto h : hits_tpc1)
            {
              if (std::find(allwire_tpc1.begin(),allwire_tpc1.end(),h) == allwire_tpc1.end()) std::cout << "PROBLEM on wire " << h << std::endl;
            }
          if (count_wires1 != 0) {
              eff1->Fill((float)count_hits1 / (float)count_wires1);
            }
        }
      Int_t count_hits3 = 0;
      Int_t count_wires3 = 0;
      Int_t halfdiff = (tpc5minwire-tpc3maxwire)/2;
      Int_t tpc3endwire = tpc3maxwire+halfdiff;
      if (tpc5minwire == 999) tpc3endwire = 112;
      Float_t thiseff3 = 1;
      if (*nwiresTPC3 != 0) {
          if (/*halfdiff>=0 &&*/ tpc3minwire != 999 && tpc3maxwire != -999)
            {
              for (int iw = 0; iw < 112; ++iw)
                {
                  if (std::find(badchannels.begin(),badchannels.end(),iw+912) != badchannels.end()) hits3string += "_";
                  else if (std::find(edgechannels.begin(),edgechannels.end(),iw+912) != edgechannels.end()) hits3string += "|";
                  else if (std::find(tpc3hits.begin(),tpc3hits.end(),iw) != tpc3hits.end()) hits3string += "X";
                  else if (std::find(currentchannels.begin(),currentchannels.end(),iw+912) == currentchannels.end()) hits3string += " ";
                  else hits3string += ".";
                }
              for (int iw = 0; iw <= tpc3endwire; ++iw)
                {
                  if ((std::find(badchannels.begin(),badchannels.end(),iw+912) != badchannels.end()) ||
                      (std::find(edgechannels.begin(),edgechannels.end(),iw+912) != edgechannels.end())) continue;
                  if (std::find(currentchannels.begin(),currentchannels.end(),iw+912) == currentchannels.end()) continue;
                  if (std::find(tpc3hits.begin(),tpc3hits.end(),iw) != tpc3hits.end()) {
                      ++count_hits3;
                      Put(tpc_wirehitsNnumerator,3,iw);
                    }
                  ++count_wires3;
                  Put(tpc_wirehitsNdenominator,3,iw);
                }
              if (count_wires3 != 0) {
                  eff3->Fill((float)count_hits3 / (float)count_wires3);
                  thiseff3 = (float)count_hits3 / (float)count_wires3;
                }
            }
        }
      Int_t count_hits5 = 0;
      Int_t count_wires5 = 0;
      Int_t tpc5beginwire = (tpc5minwire-halfdiff<0) ? 0 : tpc5minwire-halfdiff;
      if (*nwiresTPC5 != 0) {
          //if (tpc5minwire-tpc3maxwire>0 && tpc3minwire != 999 && tpc3maxwire != -999 && tpc5minwire != 999 && tpc5maxwire != -999)
          if (/*halfdiff>0 &&*/ tpc5minwire != 999 && tpc5maxwire != -999)
            {
              for (int iw = 0; iw < 112; ++iw)
                {
                  if (std::find(badchannels.begin(),badchannels.end(),iw+1424) != badchannels.end()) hits5string += "_";
                  else if (std::find(edgechannels.begin(),edgechannels.end(),iw+1424) != edgechannels.end()) hits5string += "|";
                  else if (std::find(tpc5hits.begin(),tpc5hits.end(),iw) != tpc5hits.end()) hits5string += "X";
                  else if (std::find(currentchannels.begin(),currentchannels.end(),iw+1424) == currentchannels.end()) hits5string += " ";
                  else hits5string += ".";
                }
              for (int iw = tpc5beginwire; iw < 112; ++iw)
                {
                  if ((std::find(badchannels.begin(),badchannels.end(),iw+1424) != badchannels.end()) ||
                      (std::find(edgechannels.begin(),edgechannels.end(),iw+1424) != edgechannels.end())) continue;
                  if (std::find(currentchannels.begin(),currentchannels.end(),iw+1424) == currentchannels.end()) continue;
                  if (std::find(tpc5hits.begin(),tpc5hits.end(),iw) != tpc5hits.end()) {
                      ++count_hits5;
                      Put(tpc_wirehitsNnumerator,5,iw);
                    }
                  ++count_wires5;
                  Put(tpc_wirehitsNdenominator,5,iw);
                }
              if (count_wires5 != 0) {
                  eff5->Fill((float)count_hits5 / (float)count_wires5);
                }
            }
        }
      Int_t count_hits7 = 0;
      Int_t count_wires7 = 0;
      if (*nwiresTPC7 != 0) {
          //if (tpc7minwire != 999 && tpc7maxwire != -999)
          {
            for (int iw = 0; iw < 112; ++iw)
              {
                if ((std::find(badchannels.begin(),badchannels.end(),iw+1936) != badchannels.end()) ||
                    (std::find(edgechannels.begin(),edgechannels.end(),iw+1936) != edgechannels.end())) continue;
                if (std::find(currentchannels.begin(),currentchannels.end(),iw+1936) == currentchannels.end()) continue;
                if (std::find(tpc7hits.begin(),tpc7hits.end(),iw) != tpc7hits.end()) {
                    ++count_hits7;
                    Put(tpc_wirehitsNnumerator,7,iw);
                  }
                ++count_wires7;
                Put(tpc_wirehitsNdenominator,7,iw);
              }
            if (count_wires7 != 0) {
                eff7->Fill((float)count_hits7 / (float)count_wires7);
              }
          }
        }


      //if (tpc3minwire != 999 && tpc3maxwire != -999 && tpc5minwire != 999 && tpc5maxwire != -999)
      //{
      //if (thiseff3 > 0.35 && thiseff3 < 0.45) {
      //std::cout << hits3string << std::endl;
      //std::cout << hits5string << std::endl;
      //std::cout << "halfdiff = " << halfdiff << "  tpc3endwire=" << tpc3endwire << "  tpc5beginwire=" << tpc5beginwire << std::endl;
      //std::cout << std::endl;
      //}

    }

  std::cout << "TPC1: ";
  for (auto n : nwires1) std::cout << n << " ";
  std::cout << std::endl << "TPC3: ";
  for (auto n : nwires3) std::cout << n << " ";
  std::cout << std::endl << "TPC5: ";
  for (auto n : nwires5) std::cout << n << " ";
  std::cout << std::endl << "TPC7: ";
  for (auto n : nwires7) std::cout << n << " ";
  std::cout << std::endl;

  TGaxis::SetMaxDigits(3);

  TCanvas * canvwire = new TCanvas("canvwire","",3100,1750);
  TPad * pad1 = new TPad("pad1","",0.01,0.25,0.33,0.75);
  TPad * pad5 = new TPad("pad5","",0.3333,0.51,0.66,0.99);
  TPad * pad3 = new TPad("pad3","",0.3333,0.01,0.66,0.5);
  TPad * pad7 = new TPad("pad7","",0.6666,0.25,0.99,0.75);
  pad7->Draw(); pad5->Draw(); pad3->Draw(); pad1->Draw();
  pad7->cd();
  tpc7->SetMarkerStyle(20);
  tpc7->SetMarkerSize(0.75);
  tpc7->Draw("LEP");
  tpc7->GetXaxis()->SetTitle("Wire #");
  tpc7->GetXaxis()->SetLimits(-0.5,111.5);
  tpc7->GetYaxis()->SetTitle("# Hits");
  tpc7->GetYaxis()->SetTitleOffset(1.4);
  tpc7->SetTitle("TPC 7");
  tpc7->SetStats(false);
  pad7->Update();
  for (auto b : nbox7) { b->SetY2(pad7->GetFrame()->GetY2()); b->Draw(); }
  pad5->cd();
  tpc5->SetMarkerStyle(20);
  tpc5->SetMarkerSize(0.75);
  tpc5->Draw("LEP");
  tpc5->GetXaxis()->SetTitle("Wire #");
  tpc5->GetXaxis()->SetLimits(-0.5,111.5);
  tpc5->GetYaxis()->SetTitle("# Hits");
  tpc5->GetYaxis()->SetTitleOffset(1.4);
  tpc5->SetTitle("TPC 5");
  tpc5->SetStats(false);
  pad5->Update();
  for (auto b : nbox5) { b->SetY2(pad5->GetFrame()->GetY2()); b->Draw(); }
  pad1->cd();
  tpc1->SetMarkerStyle(20);
  tpc1->SetMarkerSize(0.75);
  tpc1->Draw("LEP");
  tpc1->GetXaxis()->SetTitle("Wire #");
  tpc1->GetXaxis()->SetLimits(-0.5,111.5);
  tpc1->GetYaxis()->SetTitle("# Hits");
  tpc1->GetYaxis()->SetTitleOffset(1.4);
  tpc1->SetTitle("TPC 1");
  tpc1->SetStats(false);
  pad1->Update();
  for (auto b : nbox1) { b->SetY2(pad1->GetFrame()->GetY2()); b->Draw(); }
  pad3->cd();
  tpc3->SetMarkerStyle(20);
  tpc3->SetMarkerSize(0.75);
  tpc3->Draw("LEP");
  tpc3->GetXaxis()->SetTitle("Wire #");
  tpc3->GetXaxis()->SetLimits(-0.5,111.5);
  tpc3->GetYaxis()->SetTitle("# Hits");
  tpc3->GetYaxis()->SetTitleOffset(1.4);
  tpc3->SetTitle("TPC 3");
  tpc3->SetStats(false);
  pad3->Update();
  for (auto b : nbox3) { b->SetY2(pad3->GetFrame()->GetY2()); b->Draw(); }
  canvwire->Update();
  canvwire->SaveAs("hits_on_wires_by_tpc.png");

  TCanvas * canveff = new TCanvas("canveff","",3100,1750);
  TPad * pad1e = new TPad("pad1e","",0.01,0.25,0.33,0.75);
  TPad * pad5e = new TPad("pad5e","",0.3333,0.51,0.66,0.99);
  TPad * pad3e = new TPad("pad3e","",0.3333,0.01,0.66,0.5);
  TPad * pad7e = new TPad("pad7e","",0.6666,0.25,0.99,0.75);
  pad7e->Draw(); pad5e->Draw(); pad3e->Draw(); pad1e->Draw();
  pad7e->cd();
  eff7->SetMarkerStyle(20);
  eff7->SetMarkerSize(0.75);
  eff7->Draw("LEP");
  eff7->GetXaxis()->SetTitle("Efficiency, #epsilon_{D}");
  eff7->SetTitle("TPC 7");
  eff7->SetStats(false);
  pad5e->cd();
  eff5->SetMarkerStyle(20);
  eff5->SetMarkerSize(0.75);
  eff5->Draw("LEP");
  eff5->GetXaxis()->SetTitle("Efficiency, #epsilon_{D}");
  eff5->SetTitle("TPC 5");
  eff5->SetStats(false);
  pad1e->cd();
  eff1->SetMarkerStyle(20);
  eff1->SetMarkerSize(0.75);
  eff1->Draw("LEP");
  eff1->GetXaxis()->SetTitle("Efficiency, #epsilon_{D}");
  eff1->SetTitle("TPC 1");
  eff1->SetStats(false);
  pad3e->cd();
  eff3->SetMarkerStyle(20);
  eff3->SetMarkerSize(0.75);
  eff3->Draw("LEP");
  eff3->GetXaxis()->SetTitle("Efficiency, #epsilon_{D}");
  eff3->SetTitle("TPC 3");
  eff3->SetStats(false);
  canveff->Update();
  canveff->SaveAs("data_efficiency_by_tpc.png");


  Int_t nw1 = 0, nw3 = 0, nw5 = 0, nw7 = 0;
  for (int iw = 0; iw < 112; ++iw)
    {
      if (tpc_wirehitsNdenominator[1].find(iw) != tpc_wirehitsNdenominator[1].end())
        {
          Float_t denom = tpc_wirehitsNdenominator[1][iw];
          Float_t numer = 0;
          if (tpc_wirehitsNnumerator[1].find(iw) != tpc_wirehitsNnumerator[1].end()) numer = tpc_wirehitsNnumerator[1][iw];
          eff1wires->SetPoint(nw1,iw,numer/denom);
          eff1wires->SetPointError(nw1,0,(1/denom)*sqrt(numer*(1-(numer/denom))));
          ++nw1;
        }
      if (tpc_wirehitsNdenominator[3].find(iw) != tpc_wirehitsNdenominator[3].end())
        {
          Float_t denom = tpc_wirehitsNdenominator[3][iw];
          Float_t numer = 0;
          if (tpc_wirehitsNnumerator[3].find(iw) != tpc_wirehitsNnumerator[3].end()) numer = tpc_wirehitsNnumerator[3][iw];
          eff3wires->SetPoint(nw3,iw,numer/denom);
          eff3wires->SetPointError(nw3,0,(1/denom)*sqrt(numer*(1-(numer/denom))));
          ++nw3;
        }
      if (tpc_wirehitsNdenominator[5].find(iw) != tpc_wirehitsNdenominator[5].end())
        {
          Float_t denom = tpc_wirehitsNdenominator[5][iw];
          Float_t numer = 0;
          if (tpc_wirehitsNnumerator[5].find(iw) != tpc_wirehitsNnumerator[5].end()) numer = tpc_wirehitsNnumerator[5][iw];
          eff5wires->SetPoint(nw5,iw,numer/denom);
          eff5wires->SetPointError(nw5,0,(1/denom)*sqrt(numer*(1-(numer/denom))));
          ++nw5;
        }
      if (tpc_wirehitsNdenominator[7].find(iw) != tpc_wirehitsNdenominator[7].end())
        {
          Float_t denom = tpc_wirehitsNdenominator[7][iw];
          Float_t numer = 0;
          if (tpc_wirehitsNnumerator[7].find(iw) != tpc_wirehitsNnumerator[7].end()) numer = tpc_wirehitsNnumerator[7][iw];
          eff7wires->SetPoint(nw7,iw,numer/denom);
          eff7wires->SetPointError(nw7,0,(1/denom)*sqrt(numer*(1-(numer/denom))));
          ++nw7;
        }
    }

  std::cout << "eff7wires->N = " << eff7wires->GetN() << " eff5wires->N = " << eff5wires->GetN() << " eff3wires->N = " << eff3wires->GetN() << " eff1wires->N = " << eff1wires->GetN() << std::endl;


  TCanvas * canvwireeff = new TCanvas("canvwireeff","",3100,1750);
  TPad * pad1w = new TPad("pad1w","",0.01,0.25,0.33,0.75);
  TPad * pad5w = new TPad("pad5w","",0.3333,0.51,0.66,0.99);
  TPad * pad3w = new TPad("pad3w","",0.3333,0.01,0.66,0.5);
  TPad * pad7w = new TPad("pad7w","",0.6666,0.25,0.99,0.75);
  pad7w->Draw(); pad5w->Draw(); pad3w->Draw(); pad1w->Draw();
  pad7w->cd();
  eff7wires->SetMarkerStyle(8);
  eff7wires->SetMarkerSize(0.75);
  eff7wires->SetLineWidth(1);
  eff7wires->Draw("ap");
  eff7wires->GetXaxis()->SetTitle("Wire #");
  eff7wires->GetYaxis()->SetTitle("Efficiency, #epsilon_{D}");
  eff7wires->GetYaxis()->SetRangeUser(0,1);
  eff7wires->GetXaxis()->SetLimits(-0.5,111.5);
  eff7wires->SetTitle("TPC 7");
  for (auto b : box7) b->Draw();
  pad5w->cd();
  eff5wires->SetMarkerStyle(8);
  eff5wires->SetMarkerSize(0.75);
  eff5wires->SetLineWidth(1);
  eff5wires->Draw("ap");
  eff5wires->GetXaxis()->SetTitle("Wire #");
  eff5wires->GetYaxis()->SetTitle("Efficiency, #epsilon_{D}");
  eff5wires->GetYaxis()->SetRangeUser(0,1);
  eff5wires->GetXaxis()->SetLimits(-0.5,111.5);
  eff5wires->SetTitle("TPC 5");
  for (auto b : box5) b->Draw();
  pad1w->cd();
  eff1wires->SetMarkerStyle(8);
  eff1wires->SetMarkerSize(0.75);
  eff1wires->SetLineWidth(1);
  eff1wires->Draw("ap");
  eff1wires->GetXaxis()->SetTitle("Wire #");
  eff1wires->GetYaxis()->SetTitle("Efficiency, #epsilon_{D}");
  eff1wires->GetYaxis()->SetRangeUser(0,1);
  eff1wires->GetXaxis()->SetLimits(-0.5,111.5);
  eff1wires->SetTitle("TPC 1");
  for (auto b : box1) b->Draw();
  pad3w->cd();
  eff3wires->SetMarkerStyle(8);
  eff3wires->SetMarkerSize(0.75);
  eff3wires->SetLineWidth(1);
  eff3wires->Draw("ap");
  eff3wires->GetXaxis()->SetTitle("Wire #");
  eff3wires->GetYaxis()->SetTitle("Efficiency, #epsilon_{D}");
  eff3wires->GetYaxis()->SetRangeUser(0,1);
  eff3wires->GetXaxis()->SetLimits(-0.5,111.5);
  eff3wires->SetTitle("TPC 3");
  for (auto b : box3) b->Draw();
  canvwireeff->Update();
  canvwireeff->SaveAs("efficiency_by_wire_by_tpc.png");
}
