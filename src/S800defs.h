////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////
////////////                       GrROOT
////////////
////////////          Purpose:
////////////                   To assist in the analysis of data from
////////////                 the gretina/S800 experimental setup.
////////////                          
////////////          Current Maintainers:
////////////                 Kathrin Wimmer  (wimmer@phys.s.u-tokyo.ac.jp)
////////////                 Eric Lunderberg (lunderberg@nscl.msu.edu)
////////////
////////////          Distribution:
////////////                   Please do not redistribute this software directly.
////////////                   If someone new wants a copy of this software,
////////////                 email one of the maintainers for the download link.
////////////                   This allows us to keep track of who has the software,
////////////                 and properly distribute updates and bug fixes.
////////////                 
////////////          Suggestions:
////////////                   We view the development of the software as a collaborative
////////////                 effort, and as such, welcome and appreciate any suggestions
////////////                 for bug fixes and improvements.
////////////
////////////          Disclaimer:
////////////                 This software is provided as-is, with no warranty.
////////////                 No current or future support is guaranteed for this software.
////////////
////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
#define S800_ID 5 //s800 identifier
#define S800_CARD29_ID 8
#define S800_PHYSDATA_ID 9
#define S800_PACKET		        0x5800

#define S800_VERSION			0x0005 

#define S800_TRIGGER_PACKET		0x5801
#define S800_TOF_PACKET			0x5802
#define S800_FP_SCINT_PACKET		0x5810
#define S800_FP_IC_PACKET		0x5820
#define S800_FP_IC_ENERGY_PACKET	0x5821
#define S800_FP_IC_TIME_PACKET		0x5822
#define S800_FP_TIME_PACKET		0x5830
#define S800_FP_CRDC_PACKET		0x5840
#define S800_FP_CRDC_RAW_PACKET		0x5841
#define S800_FP_CRDC_SAMPLE_PACKET	0x5842
#define S800_FP_CRDC_PAD_PACKET		0x5843
#define S800_FP_CRDC_DSP_PACKET		0x5844
#define S800_FP_CRDC_ANODE_PACKET	0x5845
#define S800_TA_PPAC_PACKET		0x5850
#define S800_II_CRDC_PACKET		0xffff
#define S800_II_CRDC_RAW_PACKET	        0xffff
#define S800_II_CRDC_ANODE_PACKET	0xffff
#define S800_TA_PIN_PACKET		0x5860
#define S800_II_TRACK_PACKET		0x5870
#define S800_II_TRACK_RAW_PACKET	0x5871
#define S800_II_TRACK_SAMPLE_PACKET	0x5872
#define S800_II_TRACK_PAD_PACKET	0x5873
#define S800_II_TRACK_DSP_PACKET	0x5874
#define S800_OB_TRACK_PACKET		0x5880
#define S800_OB_TRACK_RAW_PACKET	0x5881
#define S800_OB_TRACK_SAMPLE_PACKET	0x5882
#define S800_OB_TRACK_PAD_PACKET	0x5883
#define S800_OB_TRACK_DSP_PACKET	0x5884
#define S800_OB_SCINT_PACKET		0x5890
#define S800_OB_PIN_PACKET		0x58A0
#define S800_FP_HODO_PACKET             0x58B0
#define S800_TIMESTAMP_PACKET		0x5803
#define S800_EVENT_NUMBER_PACKET        0x5804
#define S800_GALOTTE_PACKET             0x58D0
#define S800_LABR_PACKET                0x58E0
#define S800_LABR_DETS                  4
#define S800_MESYTDC_PACKET             0x58f0

#define S800_CRDC_MAXWIDTH              32;
#define S800_DEBUG                      0;// 0 false, 1 true

#define S800_FP_CRDC_CHANNELS           224
#define S800_FP_CRDC_SAMPLES            512
#define S800_FP_CRDC_MAX_WIDTH          32
#define S800_FP_IC_CHANNELS             16
#define S800_TRACK_PARAMETERS           6
#define S800_TRACK_COEFFICIENTS         250
#define S800_II_CRDC_CHANNELS           64
#define S800_II_TRACK_CHANNELS          256
#define S800_II_TRACK_MAX_WIDTH         32
#define S800_II_TRACK_PARAMETERS        6
#define S800_II_TRACK_COEFFICIENTS      200

#ifndef CS800_LINK_TOFTAC
#define CS800_LINK_TOFTAC
#endif
