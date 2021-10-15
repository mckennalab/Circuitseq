# PlasmidSeq with Recombinant Tn5

## Introduction

### Approach for nanopore sequencing of plasmids using recombinant Tn5 

#### Materials

* Nanopore LSK-110 kit
* NEB Quick Ligase
* Recombinant Tn5 (at 40ng/ul in 50% glycerol)
* 10x EzTn5 buffer
* 1M DTT
* NEB End Repair kit
* NEB FFPE kit
* Ampure beads
* 70% ethanol
* Nanopore Flongle
* Annealed oligo barcodes at 2.5uM

## Procedure
 
### Tagmentation reaction

This step differs depending on how many samples are being used. In general, you want 2-5ug of plasmid (total) so if you are sequencing 20 plasmids you want to put more plasmid and Tn5 per sample than if you were doing 96 plasmids
Set up the tagmentation reaction in the following order:

1. - Plasmid
2. - Barcode
3. - Master mix (Tn5, EzTn5 buffer, DTT, water)

### Tagmentation reaction 

in a PCR machine (skirted plates dont fit in our normal PCR machines, you can only use the RT if you are using a plate)

10min at 23C (this is to let Tn5 load)
10min at 37C (this is to let Tn5 tagment)
5min at 55C (this is unclear, protocols for library preps skip the 37C step and go straight to 55C, but doing 55C only doesnt work as well)
Stop the reaction by adding 1 volume of 0.2% SDS. Mix well and incubate at RT for >5 mins 
Combine all samples, add 1 vol H2O, and 0.5x ampure beads. Follow the normal ampure bead clean up. Be careful not to let the beads dry because the longer fragments of DNA will never elute. 
Elute in 25ul of H2O 


## LSK-110 library prep


This is basically the standard nanopore protocol
FFPE + End-repair
Mix well and run the following PCR protocol
7.5min at 20C
7.5min at 65C
Add 30ul of H2O, 30ul of ampure beads, follow a standard ampure clean up protocol
Elute in 30ul H2O
LSK-110 ligation, set this reaction up on ice, adding the reagents in the following order quickly, and then immediately mix well 
Allow the reaction to proceed at RT for 10min (this is a good time to do a flow cell check, see "Loading the flongle")
Add 50ul of H2O and 50ul of ampure beads
Wash 2x with 125ul ONT Long Fragment Wash Buffer 
Elute in 6-8ul of ONT Elution buffer
Loading the flongle
Thaw the FLT buffer from the LSK-110 kit and Flush buffer, Sequencing buffer, and Loading beads from a flongle-expansion kit 
Take a flongle out of the 4C
Put it into the minion + flongle adaptor and run a flow cell check
Gently peel back the sticker that is sealing the sample loading port (this is pulled back parallel to the flongle, away from the flongle ID) 
Dab away any remaining flongle storage solution (yellow stuff) from the sticker (but not from the port because you dont want to capillary action it away from the flow cell)
Prepare flush buffer (117ul of flush buffer + 3ul of flush tether)
▲
Using a P200 pipet (yellow tips work better than filter tips since they are sturdier) place the tip in the sample port and increase the volume of the pipet to gently remove 2-3ul from the flow cell (This should remove any bubbles, I have never found any in a flongle) 
Using a P200 set at 119ul add the flush buffer to the sample port without introducing bubbles
Set up the sequencing reaction, NB: shake the loadng beads well immediatly prior to loading because they settle very quickly
Load the library and sequence for 24-36 hours 
