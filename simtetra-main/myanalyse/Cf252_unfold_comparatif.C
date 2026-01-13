#ifdef __CLING__
#pragma cling optimize(0)
#endif
void Cf252_unfold_comparatif()
{
//=========Macro generated from canvas: cCompare_Unfold_Ref/Unfolded vs Reference
//=========  (Thu Dec 18 14:44:39 2025) by ROOT version 6.28/06
   TCanvas *cCompare_Unfold_Ref = new TCanvas("cCompare_Unfold_Ref", "Unfolded vs Reference",756,121,1000,800);
   gStyle->SetOptStat(0);
   cCompare_Unfold_Ref->Range(-1140.11,-26.73877,10260.99,11.60249);
   cCompare_Unfold_Ref->SetFillColor(0);
   cCompare_Unfold_Ref->SetBorderMode(0);
   cCompare_Unfold_Ref->SetBorderSize(2);
   cCompare_Unfold_Ref->SetLogy();
   cCompare_Unfold_Ref->SetGridx();
   cCompare_Unfold_Ref->SetGridy();
   cCompare_Unfold_Ref->SetFrameBorderMode(0);
   cCompare_Unfold_Ref->SetFrameBorderMode(0);
   Double_t xAxis1[156] = {0, 11, 15.28217, 20.42795, 26.47962, 33.47547, 41.45046, 50.43683, 60.46451, 71.56147, 83.75398, 97.06682, 111.5235, 127.1462, 143.9564, 161.9743, 181.2195, 201.7106, 223.4657, 246.5021, 270.8367, 296.4857, 323.4647, 351.7891, 381.4737, 412.5328, 444.9804, 478.8302, 514.0955, 550.7892, 588.9241, 628.5124, 669.5663, 712.0977, 756.118, 801.6387, 848.671, 897.2256, 947.3133, 998.9447, 1052.13, 1106.879, 1163.203, 1221.11, 1280.611, 1341.715, 1404.431, 1468.768, 1534.736, 1602.344, 1671.6, 1742.513, 1815.091, 1889.343, 1965.277, 2042.902, 2122.225, 2203.255, 2285.999, 2370.466, 2456.662, 2544.596, 2634.275, 2725.706, 2818.897, 2913.855, 3010.588, 3109.101, 3209.404, 3311.501, 3415.401, 3521.11, 3628.635, 3737.982, 3849.158, 3962.17, 4077.024, 4193.727, 4312.285, 4432.704, 4554.99, 4679.15, 4805.19, 4933.115, 5062.933, 5194.649, 5328.268, 5463.797, 5601.241, 5740.607, 5881.9, 6025.126, 6170.29, 6317.397, 6466.455, 6617.467, 6770.44, 6925.379, 7082.289, 7241.175, 7402.044, 7564.9, 7729.748, 7896.594, 8065.443, 8236.3, 8409.169, 8584.057, 8760.968, 8939.907, 9120.878, 9303.888, 9488.941, 9676.041, 9865.193, 10056.4, 10249.68, 10445.01, 10642.42, 10841.91, 11043.48, 11247.13, 11452.87, 11660.71, 11870.65, 12082.68, 12296.83, 12513.09, 12731.47, 12951.97, 13174.59, 13399.35, 13626.23, 13855.26, 14086.43, 14319.75, 14555.22, 14792.84, 15032.63, 15274.58, 15518.69, 15764.98, 16013.45, 16264.09, 16516.93, 16771.94, 17029.16, 17288.57, 17550.17, 17813.99, 18080.01, 18348.25, 18618.7, 18891.37, 19166.26, 19443.39}; 
   
   TH1D *hUnfold_Bayes_perFission_clone__1 = new TH1D("hUnfold_Bayes_perFission_clone__1","Cf-252 neutron spectrum",155, xAxis1);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(2,8.213102e-177);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(3,0.00416219);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(4,0.008424778);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(5,0.009938674);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(6,0.005150044);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(8,0.008858217);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(9,0.005531068);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(10,0.004600899);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(11,0.002841378);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(12,0.001874755);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(13,0.001683126);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(14,0.001200595);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(15,0.0007635492);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(16,0.0006629761);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(17,0.001046768);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(18,0.001383442);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(19,0.001350862);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(20,0.001683482);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(21,0.001830744);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(22,0.001635444);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(23,0.001597061);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(24,0.001540912);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(25,0.001176569);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(26,0.001321357);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(27,0.0009940218);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(28,0.00187073);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(29,0.0006746399);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(30,0.0008506588);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(31,0.0006139805);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(32,0.0007993223);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(33,0.000638829);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(34,0.0006818038);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(35,0.0006014577);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(36,0.0005572134);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(37,0.0005590587);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(38,0.0004760517);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(39,0.0005270748);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(40,0.000415827);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(41,0.0004870987);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(42,0.0007735525);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(43,0.0003759126);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(44,0.0003472192);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(45,0.0003172639);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(46,0.0003160197);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(47,0.0002810568);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(48,0.0002857411);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(49,0.0002826579);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(50,0.0002109419);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(51,0.0002549533);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(52,0.0003913506);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(53,0.00018918);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(54,0.0002200832);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(55,0.0001503544);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(56,0.0001851883);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(57,0.0001514348);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(58,0.0001465138);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(59,0.0001174564);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(60,0.0002001035);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(61,0.0001043552);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(62,6.788542e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(63,7.83211e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(64,6.230389e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(65,6.180451e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(66,6.685116e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(67,3.009071e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(68,0.0001033321);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(69,5.452776e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(70,5.077478e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(71,3.00058e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(72,3.15076e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(73,1.590954e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(74,2.620187e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(75,7.112402e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(76,3.500171e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(77,1.672694e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(78,2.168185e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(79,9.194429e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(80,1.19928e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(81,6.602008e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(82,8.654146e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(83,1.193717e-05);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(84,5.96566e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(85,6.947578e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(86,2.699233e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(87,6.148061e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(88,1.860786e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(89,6.09012e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(90,1.408679e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(91,2.030282e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(92,3.036044e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(93,3.034247e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(94,1.541223e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(95,2.367245e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(96,9.676041e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(97,3.131523e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(98,2.625423e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(99,2.390377e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(100,1.536638e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(101,6.774526e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(102,7.355101e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(103,1.138061e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(104,8.729206e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(105,7.194912e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(106,1.878467e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(107,9.52092e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(108,8.565628e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(109,9.272738e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(110,3.005442e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(111,1.060074e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(112,5.883967e-12);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(113,5.220428e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(114,1.84451e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(115,6.975374e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(116,2.308154e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(117,3.766963e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(118,5.573098e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(119,1.148225e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(120,8.931831e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(121,1.413241e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(122,1.132396e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(123,8.174585e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(124,3.744561e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(125,2.009663e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(126,2.440057e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(127,2.536887e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(128,1.571535e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(129,5.736339e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(130,5.517231e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(131,6.290636e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(132,1.954904e-13);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(133,5.816861e-20);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(134,5.901348e-23);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(135,3.974349e-26);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(136,5.352713e-26);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(137,7.088668e-22);
   hUnfold_Bayes_perFission_clone__1->SetBinContent(138,4.428096e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(2,2.064708e-93);
   hUnfold_Bayes_perFission_clone__1->SetBinError(3,1.469826e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinError(4,2.091145e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinError(5,2.271272e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinError(6,1.634972e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinError(8,2.144263e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinError(9,1.694375e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinError(10,1.545348e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinError(11,1.214422e-06);
   hUnfold_Bayes_perFission_clone__1->SetBinError(12,9.864547e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(13,9.346805e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(14,7.89411e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(15,6.295402e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(16,5.866159e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(17,7.371064e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(18,8.473941e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(19,8.373566e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(20,9.347794e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(21,9.748072e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(22,9.21346e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(23,9.1047e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(24,8.94322e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(25,7.814723e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(26,8.281614e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(27,7.182951e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(28,9.853952e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(29,5.917536e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(30,6.644811e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(31,5.645238e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(32,6.441187e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(33,5.758339e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(34,5.948872e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(35,5.587371e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(36,5.377937e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(37,5.386835e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(38,4.970865e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(39,5.230474e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(40,4.645808e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(41,5.02821e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(42,6.336506e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(43,4.417213e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(44,4.245284e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(45,4.05803e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(46,4.050065e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(47,3.81946e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(48,3.851158e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(49,3.830324e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(50,3.308921e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(51,3.63777e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(52,4.507004e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(53,3.133592e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(54,3.379858e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(55,2.793592e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(56,3.100357e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(57,2.803611e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(58,2.757682e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(59,2.469127e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(60,3.222792e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(61,2.327352e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(62,1.877125e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(63,2.01625e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(64,1.798302e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(65,1.79108e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(66,1.862771e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(67,1.249744e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(68,2.315914e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(69,1.68234e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(70,1.623413e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(71,1.24798e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(72,1.278829e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(73,9.087276e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(74,1.166195e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(75,6.075934e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(76,1.347875e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(77,9.317794e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(78,1.060848e-07);
   hUnfold_Bayes_perFission_clone__1->SetBinError(79,6.908238e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(80,7.889786e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(81,5.853868e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(82,6.702195e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(83,7.871467e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(84,5.564603e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(85,6.00512e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(86,3.743045e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(87,5.649032e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(88,3.1078e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(89,5.62235e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(90,2.704026e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(91,1.026557e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(92,1.255333e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(93,3.968536e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(94,2.82838e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(95,3.505311e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(96,2.241062e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(97,1.274919e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(98,3.691514e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(99,1.113879e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(100,2.82417e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(101,1.875186e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(102,1.953886e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(103,7.685774e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(104,6.731197e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(105,1.932492e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(106,3.122531e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(107,2.223025e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(108,2.108553e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(109,2.19386e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(110,1.24899e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(111,7.417764e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(112,5.526371e-11);
   hUnfold_Bayes_perFission_clone__1->SetBinError(113,5.205447e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(114,3.094179e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(115,6.01712e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(116,1.094554e-08);
   hUnfold_Bayes_perFission_clone__1->SetBinError(117,4.421815e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(118,5.378402e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(119,2.441285e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinError(120,6.808872e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinError(121,8.564719e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(122,7.666623e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinError(123,6.513851e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinError(124,1.394137e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(125,3.229732e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(126,3.558811e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinError(127,3.628737e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(128,9.031647e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinError(129,5.456603e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinError(130,5.351376e-09);
   hUnfold_Bayes_perFission_clone__1->SetBinError(131,5.714157e-10);
   hUnfold_Bayes_perFission_clone__1->SetBinError(132,1.007321e-11);
   hUnfold_Bayes_perFission_clone__1->SetBinError(133,5.494767e-15);
   hUnfold_Bayes_perFission_clone__1->SetBinError(134,1.750171e-16);
   hUnfold_Bayes_perFission_clone__1->SetBinError(135,4.541904e-18);
   hUnfold_Bayes_perFission_clone__1->SetBinError(136,5.270987e-18);
   hUnfold_Bayes_perFission_clone__1->SetBinError(137,6.065788e-16);
   hUnfold_Bayes_perFission_clone__1->SetBinError(138,4.79417e-09);
   hUnfold_Bayes_perFission_clone__1->SetMinimum(1.245546e-23);
   hUnfold_Bayes_perFission_clone__1->SetMaximum(5.866349e+07);
   hUnfold_Bayes_perFission_clone__1->SetEntries(155);
   hUnfold_Bayes_perFission_clone__1->SetDirectory(nullptr);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000cc");
   hUnfold_Bayes_perFission_clone__1->SetLineColor(ci);
   hUnfold_Bayes_perFission_clone__1->SetLineWidth(2);
   hUnfold_Bayes_perFission_clone__1->GetXaxis()->SetTitle("E_{true} [keV]");
   hUnfold_Bayes_perFission_clone__1->GetXaxis()->SetRange(1,110);
   hUnfold_Bayes_perFission_clone__1->GetXaxis()->SetLabelFont(42);
   hUnfold_Bayes_perFission_clone__1->GetXaxis()->SetTitleOffset(1);
   hUnfold_Bayes_perFission_clone__1->GetXaxis()->SetTitleFont(42);
   hUnfold_Bayes_perFission_clone__1->GetYaxis()->SetTitle("(arb. units or per fission)");
   hUnfold_Bayes_perFission_clone__1->GetYaxis()->SetLabelFont(42);
   hUnfold_Bayes_perFission_clone__1->GetYaxis()->SetTitleFont(42);
   hUnfold_Bayes_perFission_clone__1->GetZaxis()->SetLabelFont(42);
   hUnfold_Bayes_perFission_clone__1->GetZaxis()->SetTitleOffset(1);
   hUnfold_Bayes_perFission_clone__1->GetZaxis()->SetTitleFont(42);
   hUnfold_Bayes_perFission_clone__1->Draw("HIST");
   Double_t xAxis2[156] = {0, 11, 15.28217, 20.42795, 26.47962, 33.47547, 41.45046, 50.43683, 60.46451, 71.56147, 83.75398, 97.06682, 111.5235, 127.1462, 143.9564, 161.9743, 181.2195, 201.7106, 223.4657, 246.5021, 270.8367, 296.4857, 323.4647, 351.7891, 381.4737, 412.5328, 444.9804, 478.8302, 514.0955, 550.7892, 588.9241, 628.5124, 669.5663, 712.0977, 756.118, 801.6387, 848.671, 897.2256, 947.3133, 998.9447, 1052.13, 1106.879, 1163.203, 1221.11, 1280.611, 1341.715, 1404.431, 1468.768, 1534.736, 1602.344, 1671.6, 1742.513, 1815.091, 1889.343, 1965.277, 2042.902, 2122.225, 2203.255, 2285.999, 2370.466, 2456.662, 2544.596, 2634.275, 2725.706, 2818.897, 2913.855, 3010.588, 3109.101, 3209.404, 3311.501, 3415.401, 3521.11, 3628.635, 3737.982, 3849.158, 3962.17, 4077.024, 4193.727, 4312.285, 4432.704, 4554.99, 4679.15, 4805.19, 4933.115, 5062.933, 5194.649, 5328.268, 5463.797, 5601.241, 5740.607, 5881.9, 6025.126, 6170.29, 6317.397, 6466.455, 6617.467, 6770.44, 6925.379, 7082.289, 7241.175, 7402.044, 7564.9, 7729.748, 7896.594, 8065.443, 8236.3, 8409.169, 8584.057, 8760.968, 8939.907, 9120.878, 9303.888, 9488.941, 9676.041, 9865.193, 10056.4, 10249.68, 10445.01, 10642.42, 10841.91, 11043.48, 11247.13, 11452.87, 11660.71, 11870.65, 12082.68, 12296.83, 12513.09, 12731.47, 12951.97, 13174.59, 13399.35, 13626.23, 13855.26, 14086.43, 14319.75, 14555.22, 14792.84, 15032.63, 15274.58, 15518.69, 15764.98, 16013.45, 16264.09, 16516.93, 16771.94, 17029.16, 17288.57, 17550.17, 17813.99, 18080.01, 18348.25, 18618.7, 18891.37, 19166.26, 19443.39}; 
   
   TH1D *hUnfold_Direct_perFission__2 = new TH1D("hUnfold_Direct_perFission__2","Direct unfolded spectrum",155, xAxis2);
   hUnfold_Direct_perFission__2->SetBinContent(3,0.005416597);
   hUnfold_Direct_perFission__2->SetBinContent(4,0.0123735);
   hUnfold_Direct_perFission__2->SetBinContent(5,0.01456697);
   hUnfold_Direct_perFission__2->SetBinContent(6,0.01839984);
   hUnfold_Direct_perFission__2->SetBinContent(8,0.01061524);
   hUnfold_Direct_perFission__2->SetBinContent(9,0.008946818);
   hUnfold_Direct_perFission__2->SetBinContent(10,0.006426598);
   hUnfold_Direct_perFission__2->SetBinContent(11,0.004228899);
   hUnfold_Direct_perFission__2->SetBinContent(12,0.003009225);
   hUnfold_Direct_perFission__2->SetBinContent(13,0.002767162);
   hUnfold_Direct_perFission__2->SetBinContent(14,0.002011017);
   hUnfold_Direct_perFission__2->SetBinContent(15,0.0009544871);
   hUnfold_Direct_perFission__2->SetBinContent(16,0.0005760091);
   hUnfold_Direct_perFission__2->SetBinContent(17,0.001002865);
   hUnfold_Direct_perFission__2->SetBinContent(18,0.001652912);
   hUnfold_Direct_perFission__2->SetBinContent(19,0.001681836);
   hUnfold_Direct_perFission__2->SetBinContent(20,0.002490415);
   hUnfold_Direct_perFission__2->SetBinContent(21,0.002971915);
   hUnfold_Direct_perFission__2->SetBinContent(22,0.002812724);
   hUnfold_Direct_perFission__2->SetBinContent(23,0.002873203);
   hUnfold_Direct_perFission__2->SetBinContent(24,0.002627252);
   hUnfold_Direct_perFission__2->SetBinContent(25,0.002294752);
   hUnfold_Direct_perFission__2->SetBinContent(26,0.002111011);
   hUnfold_Direct_perFission__2->SetBinContent(27,0.002015152);
   hUnfold_Direct_perFission__2->SetBinContent(28,0.003098722);
   hUnfold_Direct_perFission__2->SetBinContent(29,0.002782154);
   hUnfold_Direct_perFission__2->SetBinContent(30,0.002494543);
   hUnfold_Direct_perFission__2->SetBinContent(31,0.002225027);
   hUnfold_Direct_perFission__2->SetBinContent(32,0.002093327);
   hUnfold_Direct_perFission__2->SetBinContent(33,0.001897761);
   hUnfold_Direct_perFission__2->SetBinContent(34,0.00166745);
   hUnfold_Direct_perFission__2->SetBinContent(35,0.001416642);
   hUnfold_Direct_perFission__2->SetBinContent(36,0.001217714);
   hUnfold_Direct_perFission__2->SetBinContent(37,0.00107713);
   hUnfold_Direct_perFission__2->SetBinContent(38,0.0009188932);
   hUnfold_Direct_perFission__2->SetBinContent(39,0.0008006468);
   hUnfold_Direct_perFission__2->SetBinContent(40,0.0007405672);
   hUnfold_Direct_perFission__2->SetBinContent(41,0.0007610448);
   hUnfold_Direct_perFission__2->SetBinContent(42,0.001414499);
   hUnfold_Direct_perFission__2->SetBinContent(43,0.001282766);
   hUnfold_Direct_perFission__2->SetBinContent(44,0.001135165);
   hUnfold_Direct_perFission__2->SetBinContent(45,0.0009652575);
   hUnfold_Direct_perFission__2->SetBinContent(46,0.0008160431);
   hUnfold_Direct_perFission__2->SetBinContent(47,0.0007210924);
   hUnfold_Direct_perFission__2->SetBinContent(48,0.0005926668);
   hUnfold_Direct_perFission__2->SetBinContent(49,0.0005167593);
   hUnfold_Direct_perFission__2->SetBinContent(50,0.0004066016);
   hUnfold_Direct_perFission__2->SetBinContent(51,0.000393371);
   hUnfold_Direct_perFission__2->SetBinContent(52,0.0007666237);
   hUnfold_Direct_perFission__2->SetBinContent(53,0.0006695667);
   hUnfold_Direct_perFission__2->SetBinContent(54,0.0005650837);
   hUnfold_Direct_perFission__2->SetBinContent(55,0.0004535796);
   hUnfold_Direct_perFission__2->SetBinContent(56,0.0003409253);
   hUnfold_Direct_perFission__2->SetBinContent(57,0.0002775439);
   hUnfold_Direct_perFission__2->SetBinContent(58,0.0001936513);
   hUnfold_Direct_perFission__2->SetBinContent(59,0.0001986451);
   hUnfold_Direct_perFission__2->SetBinContent(60,0.0003014585);
   hUnfold_Direct_perFission__2->SetBinContent(61,0.0003071155);
   hUnfold_Direct_perFission__2->SetBinContent(62,0.0002481023);
   hUnfold_Direct_perFission__2->SetBinContent(63,0.0002088008);
   hUnfold_Direct_perFission__2->SetBinContent(64,0.0001648805);
   hUnfold_Direct_perFission__2->SetBinContent(65,0.0001225835);
   hUnfold_Direct_perFission__2->SetBinContent(66,0.0001078886);
   hUnfold_Direct_perFission__2->SetBinContent(67,9.078623e-05);
   hUnfold_Direct_perFission__2->SetBinContent(68,0.0001739901);
   hUnfold_Direct_perFission__2->SetBinContent(69,0.0001418764);
   hUnfold_Direct_perFission__2->SetBinContent(70,0.0001119212);
   hUnfold_Direct_perFission__2->SetBinContent(71,6.833732e-05);
   hUnfold_Direct_perFission__2->SetBinContent(72,4.655684e-05);
   hUnfold_Direct_perFission__2->SetBinContent(73,3.764251e-05);
   hUnfold_Direct_perFission__2->SetBinContent(74,6.133513e-05);
   hUnfold_Direct_perFission__2->SetBinContent(75,5.481968e-05);
   hUnfold_Direct_perFission__2->SetBinContent(76,5.035006e-05);
   hUnfold_Direct_perFission__2->SetBinContent(77,3.638176e-05);
   hUnfold_Direct_perFission__2->SetBinContent(78,2.574582e-05);
   hUnfold_Direct_perFission__2->SetBinContent(79,1.52386e-05);
   hUnfold_Direct_perFission__2->SetBinContent(80,1.817117e-05);
   hUnfold_Direct_perFission__2->SetBinContent(81,3.91237e-05);
   hUnfold_Direct_perFission__2->SetBinContent(82,2.507557e-05);
   hUnfold_Direct_perFission__2->SetBinContent(83,2.116752e-05);
   hUnfold_Direct_perFission__2->SetBinContent(84,6.901102e-06);
   hUnfold_Direct_perFission__2->SetBinContent(85,3.888644e-06);
   hUnfold_Direct_perFission__2->SetBinContent(86,1.203224e-05);
   hUnfold_Direct_perFission__2->SetBinContent(87,1.37697e-05);
   hUnfold_Direct_perFission__2->SetBinContent(88,1.000222e-05);
   hUnfold_Direct_perFission__2->SetBinContent(89,6.177032e-06);
   hUnfold_Direct_perFission__2->SetBinContent(90,2.925267e-06);
   hUnfold_Direct_perFission__2->SetBinContent(91,4.429491e-07);
   hUnfold_Direct_perFission__2->SetBinContent(92,4.203286e-06);
   hUnfold_Direct_perFission__2->SetBinContent(93,4.632516e-06);
   hUnfold_Direct_perFission__2->SetBinContent(94,3.707707e-06);
   hUnfold_Direct_perFission__2->SetBinContent(95,3.311511e-06);
   hUnfold_Direct_perFission__2->SetBinContent(96,1.651808e-06);
   hUnfold_Direct_perFission__2->SetBinContent(97,1.499707e-06);
   hUnfold_Direct_perFission__2->SetBinContent(98,2.632409e-07);
   hUnfold_Direct_perFission__2->SetBinContent(99,1.174254e-06);
   hUnfold_Direct_perFission__2->SetBinContent(100,2.170777e-06);
   hUnfold_Direct_perFission__2->SetBinContent(101,1.307108e-06);
   hUnfold_Direct_perFission__2->SetBinContent(102,1.578797e-06);
   hUnfold_Direct_perFission__2->SetBinContent(103,2.306566e-07);
   hUnfold_Direct_perFission__2->SetBinContent(105,5.775699e-07);
   hUnfold_Direct_perFission__2->SetBinContent(106,3.600217e-07);
   hUnfold_Direct_perFission__2->SetBinContent(107,1.652991e-06);
   hUnfold_Direct_perFission__2->SetBinContent(110,2.612496e-07);
   hUnfold_Direct_perFission__2->SetBinContent(111,4.704545e-07);
   hUnfold_Direct_perFission__2->SetBinContent(112,6.26224e-07);
   hUnfold_Direct_perFission__2->SetBinContent(113,2.160658e-07);
   hUnfold_Direct_perFission__2->SetBinContent(115,2.279953e-07);
   hUnfold_Direct_perFission__2->SetBinContent(118,2.72847e-07);
   hUnfold_Direct_perFission__2->SetBinContent(119,2.593708e-07);
   hUnfold_Direct_perFission__2->SetBinContent(125,2.638258e-07);
   hUnfold_Direct_perFission__2->SetBinContent(128,5.510686e-08);
   hUnfold_Direct_perFission__2->SetBinContent(137,2.827092e-07);
   hUnfold_Direct_perFission__2->SetBinContent(138,1.181295e-07);
   hUnfold_Direct_perFission__2->SetBinError(3,1.676749e-06);
   hUnfold_Direct_perFission__2->SetBinError(4,2.53426e-06);
   hUnfold_Direct_perFission__2->SetBinError(5,2.749727e-06);
   hUnfold_Direct_perFission__2->SetBinError(6,3.090381e-06);
   hUnfold_Direct_perFission__2->SetBinError(8,2.347307e-06);
   hUnfold_Direct_perFission__2->SetBinError(9,2.15496e-06);
   hUnfold_Direct_perFission__2->SetBinError(10,1.826398e-06);
   hUnfold_Direct_perFission__2->SetBinError(11,1.481558e-06);
   hUnfold_Direct_perFission__2->SetBinError(12,1.249776e-06);
   hUnfold_Direct_perFission__2->SetBinError(13,1.198456e-06);
   hUnfold_Direct_perFission__2->SetBinError(14,1.021675e-06);
   hUnfold_Direct_perFission__2->SetBinError(15,7.03866e-07);
   hUnfold_Direct_perFission__2->SetBinError(16,5.467888e-07);
   hUnfold_Direct_perFission__2->SetBinError(17,7.214831e-07);
   hUnfold_Direct_perFission__2->SetBinError(18,9.262534e-07);
   hUnfold_Direct_perFission__2->SetBinError(19,9.343225e-07);
   hUnfold_Direct_perFission__2->SetBinError(20,1.136949e-06);
   hUnfold_Direct_perFission__2->SetBinError(21,1.242004e-06);
   hUnfold_Direct_perFission__2->SetBinError(22,1.208282e-06);
   hUnfold_Direct_perFission__2->SetBinError(23,1.221204e-06);
   hUnfold_Direct_perFission__2->SetBinError(24,1.167766e-06);
   hUnfold_Direct_perFission__2->SetBinError(25,1.091372e-06);
   hUnfold_Direct_perFission__2->SetBinError(26,1.046767e-06);
   hUnfold_Direct_perFission__2->SetBinError(27,1.022725e-06);
   hUnfold_Direct_perFission__2->SetBinError(28,1.268225e-06);
   hUnfold_Direct_perFission__2->SetBinError(29,1.201699e-06);
   hUnfold_Direct_perFission__2->SetBinError(30,1.13789e-06);
   hUnfold_Direct_perFission__2->SetBinError(31,1.074664e-06);
   hUnfold_Direct_perFission__2->SetBinError(32,1.042374e-06);
   hUnfold_Direct_perFission__2->SetBinError(33,9.92489e-07);
   hUnfold_Direct_perFission__2->SetBinError(34,9.30318e-07);
   hUnfold_Direct_perFission__2->SetBinError(35,8.575018e-07);
   hUnfold_Direct_perFission__2->SetBinError(36,7.950191e-07);
   hUnfold_Direct_perFission__2->SetBinError(37,7.477199e-07);
   hUnfold_Direct_perFission__2->SetBinError(38,6.906173e-07);
   hUnfold_Direct_perFission__2->SetBinError(39,6.446521e-07);
   hUnfold_Direct_perFission__2->SetBinError(40,6.199935e-07);
   hUnfold_Direct_perFission__2->SetBinError(41,6.285069e-07);
   hUnfold_Direct_perFission__2->SetBinError(42,8.568528e-07);
   hUnfold_Direct_perFission__2->SetBinError(43,8.159785e-07);
   hUnfold_Direct_perFission__2->SetBinError(44,7.675991e-07);
   hUnfold_Direct_perFission__2->SetBinError(45,7.07826e-07);
   hUnfold_Direct_perFission__2->SetBinError(46,6.508209e-07);
   hUnfold_Direct_perFission__2->SetBinError(47,6.117872e-07);
   hUnfold_Direct_perFission__2->SetBinError(48,5.546388e-07);
   hUnfold_Direct_perFission__2->SetBinError(49,5.179038e-07);
   hUnfold_Direct_perFission__2->SetBinError(50,4.593984e-07);
   hUnfold_Direct_perFission__2->SetBinError(51,4.518623e-07);
   hUnfold_Direct_perFission__2->SetBinError(52,6.308063e-07);
   hUnfold_Direct_perFission__2->SetBinError(53,5.895245e-07);
   hUnfold_Direct_perFission__2->SetBinError(54,5.415784e-07);
   hUnfold_Direct_perFission__2->SetBinError(55,4.852122e-07);
   hUnfold_Direct_perFission__2->SetBinError(56,4.206632e-07);
   hUnfold_Direct_perFission__2->SetBinError(57,3.795516e-07);
   hUnfold_Direct_perFission__2->SetBinError(58,3.170408e-07);
   hUnfold_Direct_perFission__2->SetBinError(59,3.211027e-07);
   hUnfold_Direct_perFission__2->SetBinError(60,3.955657e-07);
   hUnfold_Direct_perFission__2->SetBinError(61,3.9926e-07);
   hUnfold_Direct_perFission__2->SetBinError(62,3.588561e-07);
   hUnfold_Direct_perFission__2->SetBinError(63,3.292085e-07);
   hUnfold_Direct_perFission__2->SetBinError(64,2.925429e-07);
   hUnfold_Direct_perFission__2->SetBinError(65,2.522441e-07);
   hUnfold_Direct_perFission__2->SetBinError(66,2.366425e-07);
   hUnfold_Direct_perFission__2->SetBinError(67,2.170776e-07);
   hUnfold_Direct_perFission__2->SetBinError(68,3.005158e-07);
   hUnfold_Direct_perFission__2->SetBinError(69,2.713689e-07);
   hUnfold_Direct_perFission__2->SetBinError(70,2.410244e-07);
   hUnfold_Direct_perFission__2->SetBinError(71,1.883362e-07);
   hUnfold_Direct_perFission__2->SetBinError(72,1.554521e-07);
   hUnfold_Direct_perFission__2->SetBinError(73,1.397797e-07);
   hUnfold_Direct_perFission__2->SetBinError(74,1.784266e-07);
   hUnfold_Direct_perFission__2->SetBinError(75,1.686837e-07);
   hUnfold_Direct_perFission__2->SetBinError(76,1.616609e-07);
   hUnfold_Direct_perFission__2->SetBinError(77,1.37419e-07);
   hUnfold_Direct_perFission__2->SetBinError(78,1.156001e-07);
   hUnfold_Direct_perFission__2->SetBinError(79,8.893596e-08);
   hUnfold_Direct_perFission__2->SetBinError(80,9.711726e-08);
   hUnfold_Direct_perFission__2->SetBinError(81,1.425033e-07);
   hUnfold_Direct_perFission__2->SetBinError(82,1.140855e-07);
   hUnfold_Direct_perFission__2->SetBinError(83,1.04819e-07);
   hUnfold_Direct_perFission__2->SetBinError(84,5.985e-08);
   hUnfold_Direct_perFission__2->SetBinError(85,4.492665e-08);
   hUnfold_Direct_perFission__2->SetBinError(86,7.902748e-08);
   hUnfold_Direct_perFission__2->SetBinError(87,8.454096e-08);
   hUnfold_Direct_perFission__2->SetBinError(88,7.205318e-08);
   hUnfold_Direct_perFission__2->SetBinError(89,5.662326e-08);
   hUnfold_Direct_perFission__2->SetBinError(90,3.896617e-08);
   hUnfold_Direct_perFission__2->SetBinError(91,1.516288e-08);
   hUnfold_Direct_perFission__2->SetBinError(92,4.670887e-08);
   hUnfold_Direct_perFission__2->SetBinError(93,4.903581e-08);
   hUnfold_Direct_perFission__2->SetBinError(94,4.386899e-08);
   hUnfold_Direct_perFission__2->SetBinError(95,4.145892e-08);
   hUnfold_Direct_perFission__2->SetBinError(96,2.928092e-08);
   hUnfold_Direct_perFission__2->SetBinError(97,2.790025e-08);
   hUnfold_Direct_perFission__2->SetBinError(98,1.168912e-08);
   hUnfold_Direct_perFission__2->SetBinError(99,2.4688e-08);
   hUnfold_Direct_perFission__2->SetBinError(100,3.3567e-08);
   hUnfold_Direct_perFission__2->SetBinError(101,2.604718e-08);
   hUnfold_Direct_perFission__2->SetBinError(102,2.862649e-08);
   hUnfold_Direct_perFission__2->SetBinError(103,1.094178e-08);
   hUnfold_Direct_perFission__2->SetBinError(105,1.731439e-08);
   hUnfold_Direct_perFission__2->SetBinError(106,1.367002e-08);
   hUnfold_Direct_perFission__2->SetBinError(107,2.929141e-08);
   hUnfold_Direct_perFission__2->SetBinError(110,1.164482e-08);
   hUnfold_Direct_perFission__2->SetBinError(111,1.562657e-08);
   hUnfold_Direct_perFission__2->SetBinError(112,1.802892e-08);
   hUnfold_Direct_perFission__2->SetBinError(113,1.059005e-08);
   hUnfold_Direct_perFission__2->SetBinError(115,1.087847e-08);
   hUnfold_Direct_perFission__2->SetBinError(118,1.190048e-08);
   hUnfold_Direct_perFission__2->SetBinError(119,1.160287e-08);
   hUnfold_Direct_perFission__2->SetBinError(125,1.170209e-08);
   hUnfold_Direct_perFission__2->SetBinError(128,5.348202e-09);
   hUnfold_Direct_perFission__2->SetBinError(137,1.211365e-08);
   hUnfold_Direct_perFission__2->SetBinError(138,7.830402e-09);
   hUnfold_Direct_perFission__2->SetMinimum(1.245546e-23);
   hUnfold_Direct_perFission__2->SetMaximum(5.866349e+07);
   hUnfold_Direct_perFission__2->SetEntries(114);
   hUnfold_Direct_perFission__2->SetDirectory(nullptr);

   ci = TColor::GetColor("#cc0000");
   hUnfold_Direct_perFission__2->SetLineColor(ci);
   hUnfold_Direct_perFission__2->SetLineStyle(2);
   hUnfold_Direct_perFission__2->SetLineWidth(2);
   hUnfold_Direct_perFission__2->GetXaxis()->SetTitle("E_{true}");
   hUnfold_Direct_perFission__2->GetXaxis()->SetRange(0,110);
   hUnfold_Direct_perFission__2->GetXaxis()->SetLabelFont(42);
   hUnfold_Direct_perFission__2->GetXaxis()->SetTitleOffset(1);
   hUnfold_Direct_perFission__2->GetXaxis()->SetTitleFont(42);
   hUnfold_Direct_perFission__2->GetYaxis()->SetTitle("Counts");
   hUnfold_Direct_perFission__2->GetYaxis()->SetLabelFont(42);
   hUnfold_Direct_perFission__2->GetYaxis()->SetTitleFont(42);
   hUnfold_Direct_perFission__2->GetZaxis()->SetLabelFont(42);
   hUnfold_Direct_perFission__2->GetZaxis()->SetTitleOffset(1);
   hUnfold_Direct_perFission__2->GetZaxis()->SetTitleFont(42);
   hUnfold_Direct_perFission__2->Draw("HIST SAME");
   
   Double_t gRef_SEB347_fx1001[408] = { 103, 109, 114, 120, 126, 132, 138, 145, 151, 158, 164, 172, 179, 186, 193, 201, 209,
   217, 224, 233, 241, 250, 258, 267, 276, 285, 294, 304, 313, 323, 332, 343, 353,
   363, 373, 384, 395, 406, 417, 428, 439, 451, 462, 474, 486, 499, 511, 524, 536,
   549, 561, 574, 587, 601, 614, 628, 642, 656, 670, 684, 698, 713, 728, 743, 757,
   773, 788, 804, 819, 835, 851, 867, 883, 900, 916, 933, 950, 967, 984, 1002, 1019,
   1037, 1054, 1072, 1090, 1109, 1127, 1146, 1164, 1183, 1202, 1221, 1240, 1260, 1279, 1299, 1319,
   1339, 1359, 1380, 1400, 1421, 1442, 1463, 1484, 1506, 1527, 1549, 1570, 1592, 1614, 1637, 1659,
   1682, 1704, 1727, 1750, 1773, 1796, 1820, 1843, 1867, 1891, 1915, 1939, 1964, 1988, 2013, 2038,
   2088, 2139, 2190, 2242, 2295, 2348, 2402, 2456, 2511, 2567, 2623, 2680, 2738, 2796, 2855, 2914,
   2974, 3035, 3096, 3158, 3221, 3284, 3347, 3412, 3477, 3542, 3608, 3675, 3743, 3811, 3879, 3949,
   4019, 4089, 4160, 4232, 4304, 4377, 4451, 4525, 4600, 4676, 4752, 4828, 4906, 4984, 5062, 5141,
   5221, 5302, 5383, 5464, 5547, 5629, 5713, 5797, 5882, 5967, 6053, 6140, 6227, 6315, 6403, 6492,
   6582, 6672, 6763, 6855, 6947, 1855, 1871, 1887, 1903, 1919, 1935, 1952, 1968, 1985, 2001, 2034,
   2068, 2101, 2136, 2170, 2205, 2240, 2275, 2311, 2347, 2383, 2419, 2456, 2493, 2531, 2569, 2607,
   2645, 2684, 2723, 2762, 2801, 2841, 2881, 2922, 2963, 3004, 3048, 3093, 3138, 3183, 3229, 3275,
   3322, 3369, 3416, 3463, 3511, 3559, 3607, 3656, 3705, 3754, 3804, 3854, 3904, 3955, 4006, 4057,
   4109, 4161, 4213, 4266, 4319, 4372, 4426, 4480, 4534, 4588, 4643, 4699, 4754, 4810, 4866, 4923,
   4980, 5037, 5094, 5152, 5210, 5269, 5328, 5387, 5446, 5506, 5566, 5627, 5687, 5748, 5810, 5871,
   5934, 5996, 6059, 6122, 6185, 6249, 6313, 6377, 6442, 6507, 6572, 6638, 6704, 6770, 6837, 6904,
   6971, 7038, 7106, 7175, 7240, 7306, 7371, 7437, 7504, 4957, 4994, 5031, 5069, 5106, 5144, 5182,
   5220, 5258, 5296, 5334, 5373, 5411, 5450, 5489, 5528, 5567, 5606, 5646, 5685, 5725, 5765, 5804,
   5844, 5885, 5925, 5965, 6006, 6047, 6087, 6128, 6169, 6211, 6252, 6293, 6335, 6377, 6419, 6461,
   6503, 6545, 6587, 6630, 6672, 6715, 6758, 6801, 6844, 6888, 6931, 6975, 7018, 7062, 7106, 7150,
   7194, 7239, 7283, 7328, 7372, 7417, 7462, 7507, 7553, 7598, 7644, 7689, 7735, 7781, 7827, 7873,
   7919, 7966, 8012, 8059, 8106, 8153, 8200, 8247, 8294, 8342, 8389, 8437, 8485, 8533, 8581, 8629,
   8678, 8726, 8775, 8823, 8872, 8921, 8970 };
   Double_t gRef_SEB347_fy1001[408] = { 0.61236, 1.0786, 0.95743, 1.2051, 1.4854, 1.5015, 1.8341, 2.2166, 2.2606, 2.6479, 2.7973, 3.1913, 3.7813, 4.9534, 6.1443, 8.9214, 16.398,
   14.489, 14.331, 14.593, 12.397, 10.229, 8.8219, 8.89, 8.4972, 9.4235, 7.8383, 8.7514, 8.7291, 11.105, 10.997, 10.336, 9.6274,
   9.6226, 17.854, 14.746, 9.3205, 9.4769, 9.7396, 9.7581, 9.1491, 8.983, 8.4052, 9.5278, 8.8638, 9.9672, 8.6132, 8.2452, 7.1809,
   7.2502, 9.0243, 10.393, 8.8295, 8.2893, 8.3955, 7.5855, 7.5507, 6.4843, 6.7958, 6.683, 7.3458, 7.0471, 5.9992, 5.5616, 5.7337,
   5.0916, 3.9743, 4.7357, 4.1931, 3.8747, 3.8584, 2.9656, 3.7006, 3.264, 2.9481, 3.4771, 3.1852, 2.9339, 2.8314, 2.6022, 2.6632,
   2.4053, 1.8623, 2.0384, 2.0268, 2.0748, 2.075, 1.9491, 1.9317, 1.8471, 1.8233, 1.7758, 1.6942, 1.6906, 1.6467, 1.6, 1.5567,
   1.572, 1.3769, 1.4813, 1.3571, 1.3906, 1.3755, 1.2846, 1.2371, 1.2055, 1.1895, 1.1777, 1.151, 1.1025, 1.0773, 1.0778, 1.0465,
   1.0524, 1.0573, 1.0202, 0.98392, 0.97938, 0.95513, 0.85438, 0.91935, 0.87298, 0.8407, 0.84673, 0.82367, 0.81279, 0.74011, 0.73147, 0.75833,
   0.74897, 0.69592, 0.65733, 0.60902, 0.56982, 0.53858, 0.49853, 0.46005, 0.42998, 0.4056, 0.37813, 0.35022, 0.32226, 0.29823, 0.27584, 0.25662,
   0.23691, 0.22051, 0.20297, 0.18431, 0.16975, 0.15774, 0.14576, 0.13503, 0.12599, 0.11547, 0.10489, 0.095547, 0.086849, 0.080405, 0.073096, 0.066441,
   0.061825, 0.056465, 0.051475, 0.047125, 0.042668, 0.038131, 0.034849, 0.032115, 0.028927, 0.026217, 0.024141, 0.021654, 0.019539, 0.017696, 0.016257, 0.014447,
   0.012915, 0.011568, 0.010599, 0.0094351, 0.0085086, 0.0077188, 0.0069298, 0.0062526, 0.0056396, 0.0050624, 0.0045312, 0.0040786, 0.003634, 0.0032195, 0.0029152, 0.0026253,
   0.0023241, 0.0020819, 0.0018785, 0.0016774, 0.0015155, 0.86661, 0.83237, 0.84106, 0.79301, 0.80037, 0.78226, 0.72672, 0.73369, 0.73825, 0.76505, 0.73979,
   0.71805, 0.6838, 0.65706, 0.63365, 0.61356, 0.62212, 0.59481, 0.58588, 0.56365, 0.57139, 0.54461, 0.52936, 0.5013, 0.48124, 0.46809, 0.45363,
   0.42661, 0.4046, 0.38926, 0.37478, 0.3585, 0.33885, 0.31758, 0.30441, 0.29288, 0.27043, 0.24557, 0.23049, 0.2211, 0.20699, 0.19992, 0.18679,
   0.17619, 0.16854, 0.1596, 0.14908, 0.13817, 0.13126, 0.12363, 0.11624, 0.10881, 0.10365, 0.095248, 0.091134, 0.08559, 0.081195, 0.077212, 0.073825,
   0.069918, 0.06607, 0.064388, 0.058269, 0.057156, 0.052351, 0.049349, 0.047078, 0.044853, 0.041793, 0.038801, 0.036162, 0.034039, 0.031656, 0.029606, 0.027621,
   0.02568, 0.024054, 0.022478, 0.020622, 0.019143, 0.018026, 0.016752, 0.015511, 0.014295, 0.013402, 0.012388, 0.011412, 0.010622, 0.0096595, 0.0089999, 0.0083096,
   0.0076444, 0.0071385, 0.0066074, 0.0061516, 0.0056182, 0.0052484, 0.0048511, 0.0044982, 0.0040379, 0.003797, 0.003521, 0.003179, 0.002941, 0.0026844, 0.0024835, 0.0023036,
   0.0021111, 0.0019627, 0.0017587, 0.0016537, 0.0016166, 0.0014978, 0.0013721, 0.0012487, 0.0011644, 0.018491, 0.017512, 0.016503, 0.015692, 0.01499, 0.014099, 0.013533,
   0.012829, 0.012207, 0.011689, 0.011008, 0.010512, 0.0099798, 0.0093962, 0.0089857, 0.0085983, 0.0081789, 0.0077092, 0.0072839, 0.0070781, 0.006635, 0.0063594, 0.006082,
   0.0056632, 0.0053884, 0.0051771, 0.0048435, 0.0045289, 0.0043619, 0.0041161, 0.0038696, 0.0036325, 0.0034348, 0.0033278, 0.0031048, 0.0029151, 0.0027568, 0.0026204, 0.0024989,
   0.0023718, 0.0022529, 0.0021228, 0.0020087, 0.0019024, 0.0018025, 0.0017035, 0.0016183, 0.0015044, 0.0014364, 0.0013631, 0.0012851, 0.0012259, 0.0011439, 0.0010875, 0.0010202,
   0.00095572, 0.0009022, 0.00085396, 0.00080092, 0.00075724, 0.00071318, 0.00067298, 0.00062674, 0.00059451, 0.00055875, 0.00052871, 0.00049806, 0.00046547, 0.00044385, 0.00041872, 0.00039777,
   0.0003716, 0.00034947, 0.00033241, 0.00031001, 0.00029264, 0.00027749, 0.00025896, 0.00024623, 0.00022857, 0.00021771, 0.00020581, 0.00019289, 0.00018067, 0.0001706, 0.00016117, 0.00015023,
   0.00014171, 0.00013213, 0.00012509, 0.00011774, 0.00010932, 0.00010328, 9.5767e-05 };
   Double_t gRef_SEB347_fex1001[408] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0 };
   Double_t gRef_SEB347_fey1001[408] = { 9.8594, 16.288, 13.5, 15.788, 17.974, 16.668, 18.526, 20.173, 18.313, 18.802, 17.066, 16.28, 15.509, 15.366, 12.921, 9.864, 2.4528,
   2.1613, 2.1321, 2.1691, 1.8385, 1.5154, 1.3064, 1.3154, 1.2554, 1.3916, 1.1562, 1.2906, 1.2872, 1.636, 1.6197, 1.522, 1.4181,
   1.4165, 2.6277, 2.1696, 1.3717, 1.3947, 1.4329, 1.4353, 1.3458, 1.3213, 1.2362, 1.4013, 1.3038, 1.466, 1.2667, 1.2126, 1.0561,
   1.0663, 1.3273, 1.92, 1.6312, 1.5314, 1.5511, 1.4014, 1.395, 1.198, 1.2556, 1.2347, 1.3572, 1.302, 1.1084, 1.0276, 1.0595,
   0.94083, 0.73422, 0.87488, 0.77465, 0.71583, 0.71281, 0.54788, 0.68366, 0.603, 0.54465, 0.64238, 0.58846, 0.54203, 0.5231, 0.48075, 0.49202,
   0.44438, 0.34406, 0.3766, 0.37446, 0.38332, 0.38336, 0.3601, 0.35688, 0.34125, 0.33687, 0.3281, 0.31301, 0.31235, 0.30424, 0.29562, 0.28762,
   0.29044, 0.25439, 0.33662, 0.30839, 0.31601, 0.31256, 0.29193, 0.28112, 0.27395, 0.27032, 0.26763, 0.26158, 0.25054, 0.24481, 0.24494, 0.23783,
   0.23917, 0.2403, 0.23187, 0.22362, 0.22259, 0.21708, 0.19419, 0.20896, 0.19842, 0.19109, 0.19247, 0.18723, 0.18476, 0.16824, 0.16628, 0.1724,
   0.17029, 0.15824, 0.14948, 0.13852, 0.12962, 0.12254, 0.11346, 0.10473, 0.1545, 0.14574, 0.13587, 0.12584, 0.11579, 0.10716, 0.099115, 0.092207,
   0.085126, 0.079235, 0.07293, 0.066227, 0.060996, 0.056679, 0.052374, 0.04852, 0.045272, 0.041491, 0.037689, 0.034333, 0.031207, 0.028892, 0.026266, 0.023874,
   0.022216, 0.02029, 0.018497, 0.016934, 0.015332, 0.013702, 0.012523, 0.01154, 0.010395, 0.0094209, 0.008675, 0.0077812, 0.0070213, 0.006359, 0.0058418, 0.0051916,
   0.0046411, 0.0041571, 0.0038089, 0.0033906, 0.0030576, 0.0027738, 0.0024903, 0.002247, 0.0020267, 0.0018193, 0.0016284, 0.0014657, 0.001306, 0.001157, 0.0010476, 0.00094348,
   0.00083526, 0.00074819, 0.00067512, 0.00060285, 0.00054468, 0.12876, 0.12461, 0.12693, 0.12068, 0.12287, 0.12119, 0.11374, 0.11599, 0.11801, 0.15035, 0.14777,
   0.14609, 0.14185, 0.13941, 0.13767, 0.13687, 0.14279, 0.14075, 0.14335, 0.14291, 0.15042, 0.14917, 0.15131, 0.14981, 0.15873, 0.15439, 0.14962,
   0.14072, 0.13346, 0.1284, 0.12362, 0.11825, 0.11177, 0.10475, 0.10041, 0.09661, 0.089203, 0.081003, 0.07603, 0.072933, 0.068278, 0.065947, 0.061617,
   0.05812, 0.055597, 0.052646, 0.049178, 0.045578, 0.0433, 0.040782, 0.038346, 0.035893, 0.034192, 0.031421, 0.030064, 0.028235, 0.026786, 0.027443, 0.02624,
   0.024851, 0.023483, 0.022886, 0.020711, 0.020315, 0.018607, 0.017541, 0.016733, 0.015942, 0.014855, 0.013791, 0.012853, 0.012099, 0.011252, 0.010523, 0.0098177,
   0.0091278, 0.0085499, 0.0079901, 0.0073303, 0.0068044, 0.0064077, 0.0059546, 0.0055136, 0.0050815, 0.0047638, 0.0044033, 0.0040567, 0.0037756, 0.0034336, 0.0031992, 0.0029538,
   0.0027173, 0.0025375, 0.0023488, 0.0021868, 0.0019971, 0.0018657, 0.0017245, 0.001599, 0.0014354, 0.0013498, 0.0012517, 0.0011301, 0.0010455, 0.0009543, 0.00088285, 0.00081893,
   0.00075048, 0.00069776, 0.00062519, 0.00058787, 0.00057471, 0.00053246, 0.00048781, 0.00044394, 0.00041395, 0.0061474, 0.005822, 0.0054867, 0.0052168, 0.0049837, 0.0046872, 0.0044992,
   0.0042651, 0.0040585, 0.0038861, 0.0036596, 0.0034947, 0.0033179, 0.0031239, 0.0029874, 0.0028586, 0.0027192, 0.002563, 0.0024216, 0.0023532, 0.0022059, 0.0021143, 0.0020221,
   0.0018828, 0.0017914, 0.0017212, 0.0016103, 0.0015057, 0.0014502, 0.0013685, 0.0012865, 0.0012077, 0.001142, 0.0011064, 0.0010322, 0.00096916, 0.00091655, 0.00087118, 0.00083079,
   0.00078856, 0.00074902, 0.00070577, 0.00066784, 0.00063249, 0.00059927, 0.00056636, 0.00053802, 0.00050017, 0.00051397, 0.00048776, 0.00045983, 0.00043867, 0.00040932, 0.00038914, 0.00036504,
   0.00034198, 0.00032283, 0.00030557, 0.00028659, 0.00027096, 0.00025519, 0.00024081, 0.00022426, 0.00021273, 0.00019993, 0.00018918, 0.00017822, 0.00016656, 0.00015882, 0.00014983, 0.00014233,
   0.00013297, 0.00012505, 0.00011894, 0.00011093, 0.00010471, 9.9292e-05, 9.2664e-05, 8.8109e-05, 8.179e-05, 7.7903e-05, 7.3645e-05, 6.9023e-05, 6.4649e-05, 6.1044e-05, 5.7672e-05, 5.3758e-05,
   5.0708e-05, 4.7279e-05, 4.476e-05, 4.2131e-05, 3.9118e-05, 3.6955e-05, 3.4268e-05 };
   TGraphErrors *gre = new TGraphErrors(408,gRef_SEB347_fx1001,gRef_SEB347_fy1001,gRef_SEB347_fex1001,gRef_SEB347_fey1001);
   gre->SetName("gRef_SEB347");
   gre->SetTitle("Reference  SEB347 LaCl3 spectrum");
   gre->SetFillStyle(1000);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.1);
   
   TH1F *Graph_gRef_SEB3471001 = new TH1F("Graph_gRef_SEB3471001","Reference  SEB347 LaCl3 spectrum",408,0,9856.7);
   Graph_gRef_SEB3471001->SetMinimum(5.53491e-05);
   Graph_gRef_SEB3471001->SetMaximum(24.62855);
   Graph_gRef_SEB3471001->SetDirectory(nullptr);
   Graph_gRef_SEB3471001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gRef_SEB3471001->SetLineColor(ci);
   Graph_gRef_SEB3471001->GetXaxis()->SetLabelFont(42);
   Graph_gRef_SEB3471001->GetXaxis()->SetTitleOffset(1);
   Graph_gRef_SEB3471001->GetXaxis()->SetTitleFont(42);
   Graph_gRef_SEB3471001->GetYaxis()->SetLabelFont(42);
   Graph_gRef_SEB3471001->GetYaxis()->SetTitleFont(42);
   Graph_gRef_SEB3471001->GetZaxis()->SetLabelFont(42);
   Graph_gRef_SEB3471001->GetZaxis()->SetTitleOffset(1);
   Graph_gRef_SEB3471001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_gRef_SEB3471001);
   
   gre->Draw("p ");
   
   TLegend *leg = new TLegend(0.55,0.68,0.88,0.88,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("hUnfold_Bayes_perFission_clone","Bayesian unfolded (per fission)","l");

   ci = TColor::GetColor("#0000cc");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hUnfold_Direct_perFission","Direct unfolded (per fission)","l");

   ci = TColor::GetColor("#cc0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(2);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("gRef_SEB347","Reference data (LaCl_{3}:Ce)","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.1);
   entry->SetTextFont(42);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.2540281,0.9346907,0.7459719,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("Cf-252 neutron spectrum");
   pt->Draw();
   cCompare_Unfold_Ref->Modified();
   cCompare_Unfold_Ref->cd();
   cCompare_Unfold_Ref->SetSelected(cCompare_Unfold_Ref);
}
