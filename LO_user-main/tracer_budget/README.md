
#### These codes are used to calculate budget terms (except for exchange flow terms) for heat/TN/DO based on hourly history files from LiveOcean, run on perigee/apogee

- get_DO_bgc_air_sea_1.py:
  
  calculate <ins>DO production</ins> from photosynthesis, <ins>DO consumption</ins> from water column remineralization & nitrification and sediment SOD, <ins>air-sea O2 exchange</ins>, and <ins>DO\*vol</ins>

  it takes around 3.5 hr to run for 1 month.
  
- get_TNTempVol_Denitri_AirSeaHeat_1.py:
  
  calculate <ins>TN\*vol</ins>, <ins>temp\*vol</ins>, <ins>air-sea heat exchange</ins>, and <ins>sediment TN loss</ins> ( = loss of TN from detritus sinking out of the bottom water column - gain of TN from NH4 generation + loss of TN from denitrification)
