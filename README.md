# GNSS-satellite-postion
Calculate satellite position from Rinex navigation file

Test on Rinex v2.10, v3.02 with GPS navigation file

# Requirements:
    numpy
    argparse
    
# Example:
  ## 
    python satpos.py --file=rinex302.18N
    
  ## Perform time correction with --timeCor=True
    python satpos.py --file=rinex210.18N --timeCor=True
  
  ## Use Householder's iteration instead of Newton's
    python satpos.py --file=rinex210.18N --iteration=Householder
