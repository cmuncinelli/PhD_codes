/**
 * @file numericalAEEIntegration.cxx
 * * @brief Numerical integration of the Ring Observable component induced by AEE.
 *
 * This script calculates the expected value of the Ring Observable component 
 * that arises artificially due to the Azimuthal Emission Efficiency (AEE) in 
 * Lambda hyperon decays. It compares an idealized flat phase space baseline 
 * against a realistic parameterization derived from experimental data weights 
 * (accounting for jet correlation \Delta\phi and \eta_\Lambda distributions).
 *
 * The AEE-induced Ring Observable component can be written as:
 * * \mathcal{R}^{AEE}_\Lambda(\theta,\theta_\Lambda,\Delta\phi) = 
 * \frac{ -\sin\theta \cos\theta_\Lambda \cos(\Delta\phi) + \cos\theta \sin\theta_\Lambda }
 * { \sqrt{ 1 - \left( \sin\theta \sin\theta_\Lambda \cos(\Delta\phi) + \cos\theta \cos\theta_\Lambda \right)^2 } }
 *
 * The integration evaluates this observable as a function of the trigger
 * (the jet) pseudorapidity (\eta_t) and calculates the marginal distribution
 * as a function of \eta_\Lambda to assess detector-induced kinematic biases.
 * We found that both \eta_t and \eta_\Lambda distributions match with what
 * we see in experimental data, so the spurious effects really do seem to be
 * caused by AEE! Moreover, we also see that the AEE effect vanishes after full
 * dOmega^_Jet and dOmega^*_Lambda integration, which leads to a pure vortex-
 * polarization observable!!! Will further verify if this is correct using
 * pp data, MC, and an approach similar to event mixing in following codes.
 *
 * --> Save this file as: numericalAEEIntegration.cxx
 * --> Run in terminal with: root -l numericalAEEIntegration.cxx
 */

#include <iostream>
#include <cmath>
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"

// 1. Parameterization for the Jet correlation Delta_phi distribution
double weight_Delta_phi(double delta_phi) {
    double base = 1.0;
    double near_side = 2.0 * exp(-0.5 * pow(delta_phi / 0.4, 2));
    double dphi_away = fabs(delta_phi) - M_PI;
    double away_side = 1.0 * exp(-0.5 * pow(dphi_away / 0.6, 2));
    return base + near_side + away_side;
}

// 2. Parameterization for the eta_Lambda distribution
double weight_eta_Lambda(double eta_L) {
    double base = 24.0;
    double left_peak  = 8.0 * exp(-0.5 * pow((eta_L + 0.7)/0.3, 2)); 
    double right_peak = 6.0 * exp(-0.5 * pow((eta_L - 0.7)/0.3, 2)); 
    return base + left_peak + right_peak;
}

// Actual array of counts estimated from quick simulations:
// (this is introduced by hand for this ToyModel test, and pretty "phenomenological". Do NOT take these at face-value)
const double data_2D[1764] = {0, 433, 463, 493, 467, 492, 465, 522, 495, 454, 440, 431, 411, 440, 430, 423, 472, 461, 492, 405, 477, 501, 473, 457, 450, 459, 453, 419, 450, 448, 487, 453, 415, 411, 470, 459, 450, 495, 447, 476, 493, 0, 0, 3852, 3888, 3810, 3788, 3794, 3797, 3649, 3395, 3577, 3382, 3386, 3299, 3438, 3551, 3558, 3418, 3575, 3714, 3611, 3753, 3703, 3791, 3650, 3684, 3615, 3449, 3560, 3463, 3345, 3369, 3517, 3413, 3493, 3452, 3613, 3708, 3722, 3826, 3811, 3870, 0, 0, 20301, 20363, 19961, 19849, 19722, 19439, 19047, 18601, 18220, 18420, 17976, 17670, 18389, 18148, 18680, 19019, 18814, 19543, 19660, 20244, 19945, 20057, 19324, 19261, 18962, 18589, 18522, 18127, 18063, 18381, 18323, 18433, 18554, 19054, 19233, 19522, 19983, 20150, 20281, 20447, 0, 0, 62095, 61268, 61095, 60161, 59188, 58445, 57461, 55721, 54608, 54429, 54432, 53794, 54283, 54762, 55097, 56571, 57335, 58992, 59731, 60865, 60393, 60148, 58883, 57938, 57013, 55836, 54695, 54271, 53850, 54250, 54664, 55372, 55889, 57200, 57792, 58943, 59636, 60891, 61844, 61548, 0, 0, 125910, 125409, 124374, 122942, 119990, 117974, 115148, 113376, 111494, 109678, 108677, 108466, 109514, 110433, 111360, 114755, 117319, 118829, 121995, 123580, 123685, 122377, 119846, 117268, 114302, 112709, 110240, 108238, 108339, 109100, 109864, 111414, 113360, 115077, 117802, 119766, 122018, 123592, 124990, 125522, 0, 0, 199998, 199465, 197601, 193274, 190149, 187135, 182991, 178732, 174890, 172819, 171510, 171147, 171854, 173278, 176900, 181356, 185774, 190652, 194259, 196337, 196461, 194027, 190574, 185730, 180243, 176831, 173710, 171906, 170799, 171841, 173243, 176663, 178559, 182822, 186308, 190656, 193742, 196511, 198318, 199796, 0, 0, 278636, 276659, 273665, 269871, 263575, 259152, 253330, 248033, 243885, 240480, 237414, 236776, 236529, 240602, 245211, 251729, 259033, 265241, 270557, 272879, 273605, 271647, 264938, 258300, 251301, 245537, 241868, 237965, 236823, 237678, 240272, 242794, 247348, 251931, 258217, 263623, 268859, 272788, 276449, 276798, 0, 0, 357409, 355567, 351103, 345861, 338713, 331027, 324559, 317457, 311909, 308117, 304415, 303328, 305049, 308763, 315871, 321702, 331954, 340765, 349264, 354370, 353047, 348876, 342863, 332003, 322615, 315543, 309244, 305700, 303454, 304102, 306882, 311522, 317456, 325037, 331796, 338379, 344949, 350348, 353534, 357028, 0, 0, 434602, 429435, 425055, 419476, 411964, 403382, 394196, 386650, 378406, 372340, 369006, 367984, 369500, 375949, 383239, 392531, 405812, 415712, 424560, 431291, 432071, 427345, 417324, 406100, 394103, 383936, 375911, 368753, 368192, 369213, 371593, 379201, 386199, 394191, 402134, 411572, 418895, 425167, 431018, 432546, 0, 0, 506693, 503760, 497698, 490366, 481447, 471813, 460332, 450553, 441593, 434752, 431211, 430449, 433547, 440405, 449760, 461537, 474844, 489292, 501669, 506955, 508849, 502044, 490084, 476260, 461860, 448765, 439666, 432846, 431948, 432815, 437137, 443712, 451126, 459914, 469368, 481728, 489094, 496230, 502932, 506958, 0, 0, 578865, 572145, 566720, 558049, 547647, 537149, 526164, 513825, 503870, 497743, 492492, 490217, 493886, 500035, 512130, 527852, 543217, 562478, 575857, 585234, 583719, 576328, 561339, 544969, 526595, 513751, 500279, 493823, 490165, 492628, 496690, 504035, 513350, 524268, 537156, 546630, 558040, 565633, 573112, 574670, 0, 0, 643034, 641208, 632894, 623638, 610478, 597463, 585053, 572913, 562304, 554038, 548596, 547346, 550589, 558925, 572667, 590510, 610484, 632531, 647732, 659809, 660450, 651006, 631688, 612634, 590027, 574267, 559569, 551032, 546532, 549163, 555471, 563230, 573833, 585781, 598408, 610556, 621175, 631271, 639233, 642949, 0, 0, 703682, 700907, 693144, 683272, 670481, 654142, 641386, 628692, 615222, 608146, 601950, 602072, 604200, 615633, 630329, 649416, 674017, 700140, 722931, 732526, 734447, 721018, 700673, 674349, 649868, 629451, 613838, 603061, 599515, 598949, 606399, 617845, 628592, 642684, 655804, 669251, 681351, 690934, 699944, 704814, 0, 0, 751621, 746665, 738702, 728222, 714013, 697290, 684282, 669385, 658450, 648328, 641671, 641674, 645382, 655923, 671953, 698136, 725540, 757880, 781846, 792630, 793810, 778931, 756967, 726160, 697897, 673777, 657449, 648122, 642689, 642136, 646846, 659168, 671762, 684281, 698709, 713493, 726453, 737511, 747425, 751025, 0, 0, 773358, 769732, 761139, 749718, 735230, 719323, 703824, 689966, 676126, 667228, 660355, 658940, 665059, 676116, 692357, 722907, 755297, 790014, 816324, 826097, 828630, 816774, 789361, 754853, 720867, 694537, 675168, 664398, 659935, 659066, 665924, 677077, 688909, 705429, 720487, 733761, 748013, 759300, 767623, 773689, 0, 0, 768086, 765349, 757247, 745708, 731027, 714852, 700930, 685720, 673450, 663320, 657957, 658117, 663167, 674689, 692139, 722167, 759372, 795084, 820413, 818834, 820916, 818996, 792097, 757628, 722405, 691987, 674854, 660891, 655746, 657509, 663308, 674438, 687289, 701053, 716693, 732517, 745025, 755792, 765757, 768483, 0, 0, 755042, 750943, 744305, 730568, 716522, 704333, 687939, 673186, 662553, 651371, 646492, 646844, 649484, 661801, 679964, 713503, 750133, 783768, 781455, 767284, 767239, 779999, 782827, 751730, 712019, 681713, 662191, 651113, 646063, 645714, 653702, 661677, 673906, 687491, 701578, 717474, 729932, 742550, 750472, 756457, 0, 0, 737607, 732617, 724635, 711264, 698553, 686188, 670921, 658200, 646432, 636771, 630098, 629732, 634646, 645728, 666756, 700709, 736127, 755827, 786870, 1063109, 1062191, 786351, 754720, 735372, 698522, 667373, 645523, 634367, 630413, 631954, 634956, 645693, 657931, 671167, 686470, 702034, 714235, 724493, 731634, 737278, 0, 0, 719791, 714817, 705619, 696178, 684419, 669722, 655715, 642682, 631859, 620839, 614392, 614765, 618874, 631837, 652637, 683944, 720412, 730121, 1079778, 2058386, 2058282, 1077862, 730247, 721007, 685155, 651978, 631177, 618309, 615857, 617220, 623917, 632341, 641667, 654905, 671121, 684666, 695341, 706471, 713980, 719345, 0, 0, 705271, 700158, 693881, 683473, 670189, 658090, 643880, 630654, 621462, 613077, 603854, 602988, 609591, 619289, 641667, 674035, 707841, 749190, 1593575, 3417394, 3419497, 1592560, 748149, 706348, 673457, 640294, 620712, 607109, 605884, 604420, 612183, 620576, 631981, 644532, 658664, 672070, 682761, 695553, 699701, 706116, 0, 0, 696398, 693767, 685441, 675116, 664689, 651189, 636681, 625014, 613107, 604830, 598423, 598674, 601717, 611564, 634902, 667885, 695658, 777898, 1987356, 4725561, 4730045, 1985520, 776401, 694459, 666894, 634284, 614545, 600841, 596342, 598830, 603328, 612914, 623626, 638591, 650269, 662140, 676430, 683506, 691427, 696358, 0, 0, 694635, 691901, 681955, 673660, 661925, 648751, 636400, 622459, 612246, 602143, 596489, 598333, 601247, 610812, 632747, 667159, 694757, 778628, 1983437, 4716721, 4714829, 1983234, 776758, 693322, 667060, 633421, 611494, 599678, 596961, 596536, 603647, 610416, 623919, 636254, 649419, 663343, 674076, 682671, 689584, 695819, 0, 0, 699412, 697006, 689114, 677809, 666694, 653449, 641095, 628270, 616872, 608827, 599006, 600081, 606401, 614162, 636609, 670224, 701678, 743419, 1583427, 3393012, 3395988, 1579184, 743903, 703822, 669476, 637071, 616684, 604660, 600120, 600803, 605860, 615093, 629083, 641430, 652903, 667123, 679713, 689740, 696197, 700256, 0, 0, 712553, 709082, 699728, 686853, 676905, 663665, 649408, 636579, 624192, 617071, 610854, 609658, 613088, 625858, 646771, 679010, 714159, 722518, 1066375, 2035119, 2032233, 1063913, 721763, 714423, 677052, 645593, 624808, 613359, 611223, 609153, 614845, 624529, 635334, 649578, 663322, 677629, 690183, 698649, 705254, 710371, 0, 0, 725980, 721360, 712770, 702645, 691474, 675847, 664134, 648649, 637668, 628219, 622205, 622587, 626255, 636582, 657172, 689112, 725689, 745900, 774534, 1043920, 1041877, 774963, 745654, 724923, 688870, 658374, 636968, 625348, 621843, 622140, 628960, 637155, 650147, 662812, 676830, 689112, 702916, 715208, 721348, 726259, 0, 0, 740404, 736386, 728907, 716920, 707087, 690403, 675523, 661750, 651837, 640311, 634542, 633094, 637972, 649445, 668649, 699276, 738544, 773099, 767982, 752995, 751316, 766727, 771253, 738508, 700235, 670139, 651293, 637610, 635474, 635802, 640590, 650303, 663866, 675551, 690938, 705077, 717825, 729807, 736486, 739807, 0, 0, 753278, 748945, 741369, 727926, 716496, 700082, 687932, 672709, 660591, 651342, 643484, 641570, 648709, 659523, 677832, 704705, 742225, 777056, 800184, 800549, 802937, 799218, 778438, 741981, 707637, 677902, 660337, 648347, 643419, 644835, 651284, 660001, 672473, 684958, 699957, 715616, 727873, 739384, 748408, 753315, 0, 0, 752376, 749127, 741895, 727476, 716956, 700390, 687432, 673856, 658958, 648932, 644157, 641877, 648400, 659365, 674481, 701610, 735346, 770392, 792327, 805645, 806445, 793040, 767908, 735297, 703626, 675148, 659823, 648602, 642147, 644849, 649925, 658947, 670964, 686691, 701024, 714090, 726515, 737746, 747796, 751455, 0, 0, 727023, 725695, 716152, 706901, 693818, 676577, 662963, 650640, 638631, 628451, 623594, 621468, 625000, 635422, 650789, 674366, 704312, 734559, 758780, 768830, 768855, 756789, 735426, 704436, 676922, 653413, 634715, 626317, 623390, 623717, 629986, 637581, 649817, 664051, 677576, 691417, 703901, 715781, 723224, 728654, 0, 0, 679612, 677498, 668002, 658461, 647311, 633258, 619120, 606922, 596348, 587242, 582856, 579943, 583773, 593218, 606171, 626250, 653600, 677486, 697300, 707303, 710420, 697356, 677073, 652541, 627535, 607448, 594756, 584541, 580573, 581255, 586600, 596263, 606933, 619254, 632132, 645559, 655891, 667288, 675025, 680459, 0, 0, 621040, 617336, 611728, 599566, 589827, 577723, 565712, 554094, 544152, 536392, 531534, 529573, 531440, 539405, 553566, 569910, 590649, 611698, 627756, 636733, 636498, 629483, 611830, 589063, 570148, 553510, 541302, 532872, 529407, 530474, 535572, 544637, 554510, 564719, 577557, 589504, 600510, 608053, 616151, 620846, 0, 0, 556111, 553290, 546994, 539885, 529063, 518728, 506033, 497703, 485783, 479525, 474728, 473086, 476641, 483038, 492628, 509790, 524728, 542000, 556057, 564904, 563813, 556427, 541236, 526046, 509038, 494646, 484965, 476598, 473757, 475090, 480578, 487351, 496204, 505981, 518125, 527861, 537378, 546224, 552724, 555854, 0, 0, 487742, 485789, 481091, 473485, 463505, 455633, 444920, 435807, 427638, 421137, 416240, 415696, 418650, 423853, 431984, 442888, 458647, 473842, 484515, 490899, 491146, 484081, 474118, 458569, 445162, 432028, 424148, 418989, 415810, 416738, 421291, 426490, 434472, 444509, 453424, 463151, 471927, 480291, 483337, 487794, 0, 0, 415348, 414161, 407909, 401929, 395123, 385782, 379634, 371820, 364051, 358997, 356119, 354114, 354786, 360222, 368861, 378422, 389278, 400624, 408606, 415327, 415436, 409624, 400289, 388615, 378161, 367631, 361472, 356823, 354362, 355581, 359580, 363358, 371575, 379114, 387029, 395002, 401684, 409190, 414703, 414831, 0, 0, 339528, 337084, 334644, 328490, 323892, 316735, 310221, 304169, 298553, 294058, 289754, 289895, 291377, 294343, 299082, 308940, 316795, 325784, 332712, 337705, 337252, 334032, 326123, 318275, 309514, 301184, 293993, 290956, 289690, 291779, 293536, 297210, 303192, 309578, 316127, 322670, 329739, 335042, 339819, 339404, 0, 0, 260705, 261309, 256830, 255087, 249862, 244434, 239156, 234900, 229580, 226656, 224040, 223234, 225276, 227169, 232976, 237697, 243092, 248775, 254044, 259212, 257624, 256008, 250412, 242774, 237367, 230980, 227040, 225149, 223628, 225165, 226114, 229728, 233879, 239145, 244498, 249976, 253381, 256792, 260166, 261398, 0, 0, 185576, 183505, 182166, 179544, 175365, 173229, 169349, 165195, 163183, 161027, 158482, 158244, 158880, 160271, 163466, 167811, 172012, 175861, 179072, 181305, 181440, 179767, 176021, 171434, 167110, 164066, 160495, 159517, 159121, 158482, 160001, 162143, 165429, 168082, 172337, 176055, 179421, 181967, 183690, 184374, 0, 0, 112961, 112572, 111977, 110193, 108406, 105642, 104518, 102619, 100718, 98898, 97891, 97672, 97907, 99607, 100973, 103205, 105921, 108074, 110240, 110554, 111333, 109937, 108557, 106141, 103154, 100971, 99405, 98451, 98395, 97932, 98378, 100344, 101863, 103854, 105914, 108582, 110639, 111862, 113391, 113488, 0, 0, 53784, 53490, 52931, 52585, 51636, 50614, 49430, 48880, 48054, 47025, 46850, 47114, 46370, 47479, 47863, 48763, 50040, 50889, 51914, 52559, 52334, 51937, 51122, 49924, 48879, 48328, 47642, 47199, 47125, 46999, 47175, 47898, 48588, 49328, 50043, 51448, 51985, 52807, 53471, 53631, 0, 0, 16565, 16606, 16856, 16227, 15707, 15755, 15423, 15145, 15103, 14877, 14831, 14587, 14791, 14660, 15249, 15209, 15654, 15766, 16082, 16064, 16205, 15987, 16132, 15633, 15395, 15216, 14758, 14659, 14604, 14828, 14944, 14938, 15245, 15345, 15456, 16076, 16254, 16458, 16588, 16848, 0, 0, 2896, 2972, 2920, 2825, 2799, 2818, 2825, 2723, 2679, 2469, 2595, 2553, 2651, 2673, 2676, 2592, 2724, 2582, 2849, 2828, 2819, 2875, 2729, 2683, 2738, 2584, 2571, 2617, 2605, 2635, 2635, 2611, 2763, 2724, 2797, 2797, 2874, 2946, 2811, 2990, 0, 0, 341, 344, 371, 388, 338, 345, 322, 334, 327, 315, 340, 307, 325, 317, 330, 375, 321, 352, 364, 344, 355, 312, 316, 358, 312, 301, 367, 315, 284, 312, 351, 346, 329, 330, 377, 347, 353, 350, 391, 366, 0};
double weight_from_data(double delta_eta, double dphi) {
    double min_dphi = -M_PI, max_dphi = M_PI;
    double min_deta = -1.8,  max_deta =  1.8;

    if (dphi      < min_dphi || dphi      >= max_dphi) return 0.0;
    if (delta_eta < min_deta || delta_eta >= max_deta)  return 0.0;

    // 40 physical bins in each axis; offset by +1 to skip underflow slot
    int bin_x = 1 + (int)(40.0 * (dphi      - min_dphi) / (max_dphi - min_dphi));
    int bin_y = 1 + (int)(40.0 * (delta_eta - min_deta) / (max_deta - min_deta));

    // Clamp to [1, 40] (boundary values with dphi exactly = max would overshoot)
    if (bin_x > 40) bin_x = 40;
    if (bin_y > 40) bin_y = 40;

    // ROOT's stride is (N_x + 2) = 42
    int idx = bin_y * 42 + bin_x;

    return data_2D[idx];
}

void numericalAEEIntegration() {
    int n_grid = 800; 
    
    double eta_L_min = -0.9, eta_L_max = 0.9;
    double dphi_min  = -M_PI, dphi_max = M_PI;
    double eta_t_min = -0.5, eta_t_max = 0.5;
    
    double d_eta_L = (eta_L_max - eta_L_min) / n_grid;
    double d_dphi  = (dphi_max - dphi_min) / n_grid;
    double d_eta_t = (eta_t_max - eta_t_min) / n_grid;
    
    int n_scan_points = 150;

    TGraph* gr_flat  = new TGraph(); 
    TGraph* gr_param = new TGraph(); 
    TGraph* gr_flat_div  = new TGraph();
    TGraph* gr_param_div = new TGraph();
    TGraph* gr_marginal_flat  = new TGraph(n_scan_points);
    TGraph* gr_marginal_param = new TGraph(n_scan_points);
    
    // =========================================================================
    // LOOP 1: Evaluate R as a function of eta_t
    // =========================================================================
    int div_point_idx = 0;
    for (int i = 0; i < n_scan_points; i++) {
        double eta_t = eta_t_min + i * ((eta_t_max - eta_t_min) / (n_scan_points - 1));
        
        double integral_flat = 0.0, weight_sum_flat = 0.0;
        double integral_param = 0.0, weight_sum_param = 0.0;
        
        for (int j = 0; j < n_grid; j++) {
            double eta_L = eta_L_min + (j + 0.5) * d_eta_L;
            for (int m = 0; m < n_grid; m++) {
                double dphi = dphi_min + (m + 0.5) * d_dphi;
                
                double num = - sinh(eta_L) * cos(dphi) + sinh(eta_t);
                double den_sq = pow(sin(dphi), 2) + pow(sinh(eta_t), 2) + pow(sinh(eta_L), 2) 
                                - 2.0*sinh(eta_t)*sinh(eta_L)*cos(dphi);
                double den = sqrt(den_sq + 1e-6); 
                double R_kinematic = num / den;
                
                double w_flat = 1.0; 
                integral_flat += R_kinematic * w_flat;
                weight_sum_flat += w_flat;
                
                // --- Apply the REAL Data Correlation ---
                double delta_eta = eta_L - eta_t; 
                
                // Fetch the measured pair density for this specific angular separation
                    // delta_eta doesn't really appear in the equation, but it was the weight that I had in hand
                double w_data_2D = weight_from_data(delta_eta, dphi);
                
                // Calculating the total weight for the current integration value:
                    // Might be over-correcting weights (two corrections), but the overall shape is comprised in this alone
                    // If you treat both the eta and the phi distribution separately, the shape of the marginal distribution
                    // in eta_Lambda is not properly conveyed! We lack the falls in the edges otherwise.
                double w_param = weight_eta_Lambda(eta_L) * w_data_2D;
                
                integral_param += R_kinematic * w_param;
                weight_sum_param += w_param;
            }
        }
        
        double y_flat = integral_flat / weight_sum_flat;
        double y_param = integral_param / weight_sum_param;
        
        gr_flat->SetPoint(i, eta_t, y_flat);
        gr_param->SetPoint(i, eta_t, y_param);
        
        if (fabs(eta_t) > 1e-4) {
            gr_flat_div->SetPoint(div_point_idx, eta_t, y_flat / eta_t);
            gr_param_div->SetPoint(div_point_idx, eta_t, y_param / eta_t);
            div_point_idx++;
        }
    }

    // =========================================================================
    // LOOP 2: Evaluate Marginal R as a function of eta_Lambda
    // =========================================================================
    for (int i = 0; i < n_scan_points; i++) {
        double eta_L = eta_L_min + i * ((eta_L_max - eta_L_min) / (n_scan_points - 1));
        
        double integral_flat = 0.0, weight_sum_flat = 0.0;
        double integral_param = 0.0, weight_sum_param = 0.0;
        
        for (int j = 0; j < n_grid; j++) {
            double eta_t = eta_t_min + (j + 0.5) * d_eta_t;
            for (int m = 0; m < n_grid; m++) {
                double dphi = dphi_min + (m + 0.5) * d_dphi;
                
                double num = - sinh(eta_L) * cos(dphi) + sinh(eta_t);
                double den_sq = pow(sin(dphi), 2) + pow(sinh(eta_t), 2) + pow(sinh(eta_L), 2) 
                                - 2.0*sinh(eta_t)*sinh(eta_L)*cos(dphi);
                double den = sqrt(den_sq + 1e-6); 
                double R_kinematic = num / den;
                
                double w_flat = 1.0; 
                integral_flat += R_kinematic * w_flat;
                weight_sum_flat += w_flat;
                
                // --- Apply the REAL Data Correlation ---
                double delta_eta = eta_L - eta_t;
                
                double w_data_2D = weight_from_data(delta_eta, dphi);
                
                double w_param = weight_eta_Lambda(eta_L) * w_data_2D;
                
                integral_param += R_kinematic * w_param;
                weight_sum_param += w_param;
            }
        }
        
        gr_marginal_flat->SetPoint(i, eta_L, integral_flat / weight_sum_flat);
        gr_marginal_param->SetPoint(i, eta_L, integral_param / weight_sum_param);
    }
    
    // =========================================================================
    // PLOTTING
    // =========================================================================
    
    auto style_graph = [](TGraph* g_flat, TGraph* g_param, const char* title) {
        g_flat->SetTitle(title);
        g_flat->SetLineColor(kRed);
        g_flat->SetLineStyle(2); 
        g_flat->SetLineWidth(3);
        
        g_param->SetLineColor(kBlue);
        g_param->SetLineWidth(3);
    };

    // A robust function to auto-scale axes without crashing ROOT
    auto adjust_axes = [](TCanvas* c, TGraph* g1, TGraph* g2) {
        c->Update(); // Force ROOT to build the TAxis first
        double min1 = TMath::MinElement(g1->GetN(), g1->GetY());
        double max1 = TMath::MaxElement(g1->GetN(), g1->GetY());
        double min2 = TMath::MinElement(g2->GetN(), g2->GetY());
        double max2 = TMath::MaxElement(g2->GetN(), g2->GetY());
        
        double ymin = std::min(min1, min2);
        double ymax = std::max(max1, max2);
        
        // Use absolute difference for the margin!
        double margin = 0.2 * fabs(ymax - ymin);
        if (margin < 1e-5) margin = 0.1 * fabs(ymin);
        if (margin < 1e-5) margin = 0.1; // Fallback if everything is exactly 0
        
        g1->GetYaxis()->SetRangeUser(ymin - margin, ymax + margin);
        c->Modified();
        c->Update();
    };

    style_graph(gr_flat, gr_param, "Kinematic Bias vs #eta_{Jet}; #eta_{Jet}; <#it{R}>");
    style_graph(gr_flat_div, gr_param_div, "Ratio Check: <#it{R}> / #eta_{Jet}; #eta_{Jet}; <#it{R}> / #eta_{Jet}");
    style_graph(gr_marginal_flat, gr_marginal_param, "Marginal Bias vs #eta_{#Lambda}; #eta_{#Lambda}; <#it{R}>");

    // Canvas 1: R vs eta_t
    TCanvas* c1 = new TCanvas("c1", "R vs eta_t", 900, 600);
    c1->SetGrid();
    gr_flat->Draw("AL"); 
    gr_param->Draw("L SAME"); 
    adjust_axes(c1, gr_flat, gr_param);
    TLegend* leg1 = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg1->AddEntry(gr_flat, "Flat Phase Space", "l");
    leg1->AddEntry(gr_param, "Realistic Weights", "l");
    leg1->Draw();

    // Canvas 2: R/eta_t vs eta_t -- Just for debugging.
        // Don't really trust the output of this particular numerical integration, as it is pretty unstable.
        // The main conclusions can already be drafted from Canvas 1 and Canvas 2. This is just a closure of sorts.
    TCanvas* c2 = new TCanvas("c2", "R/eta_t vs eta_t", 900, 600);
    c2->SetGrid();
    gr_flat_div->Draw("AL"); 
    gr_param_div->Draw("L SAME"); 
    adjust_axes(c2, gr_flat_div, gr_param_div);
    TLegend* leg2 = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg2->AddEntry(gr_flat, "Flat Phase Space", "l");
    leg2->AddEntry(gr_param, "Realistic Weights", "l");
    leg2->Draw();

    // Canvas 3: Marginal R vs eta_Lambda
    TCanvas* c3 = new TCanvas("c3", "Marginal R vs eta_L", 900, 600);
    c3->SetGrid();
    gr_marginal_flat->Draw("AL"); 
    gr_marginal_param->Draw("L SAME"); 
    adjust_axes(c3, gr_marginal_flat, gr_marginal_param);
    TLegend* leg3 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg3->AddEntry(gr_flat, "Flat Phase Space", "l");
    leg3->AddEntry(gr_param, "Realistic Weights", "l");
    leg3->Draw();
}