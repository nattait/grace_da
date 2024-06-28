SUPPLEMENTARY DATA FOR "GLOBAL-SCALE REGIONALIZATION OF HYDROLOGIC MODEL
PARAMETERS" BY BECK ET AL.

The directories represent the ensemble members (ten in total). Each directory
contains 14 GeoTiff maps with regionalized parameters (named Par_X.tif, where X
denotes the parameter) and a GeoTiff map with the identifier of the donor
catchment (ID1.tif). The grids have a spatial resolution of 0.5°.

The Matlab implementation of HBV is also included (HBV_lumped.m). Only the 
lumped version is provided, the gridded version is available on request.

Also included is a shapefile called Catchments with metadata for all catchments
with area <10,000 km² and daily observed streamflow data. The catchment latitude
and longitude coordinates represent the centroids of the catchments rather than
the gauging station locations. The shapefile contains the following fields:
  - ID1:        Identifier assigned in our study (relates to the ID1.tif grids)
  - ID2:        Official catchment identifier from data provider
  - RivName:  	River name
  - StatName:  	Station name
  - StatLat:	Station latitude [°]
  - StatLon:	Station longitude [°]
  - Class:      "0" if catchment did not satisfy the inclusion criteria and 
                hence was excluded, "1" if catchment used as donor, "2" if used
                for evaluation
  - Area:       Calculated catchment area [km²] from Lehner (2012, Report 41
                in the GRDC Report Series)
  - RecLength:  Length of observed streamflow record during 1979–2005 [years]
  - ForGain:    Catchment fraction of forest gain [-]
  - ForLoss:    Catchment fraction of forest loss [-]
  - IrrArea:    Catchment fraction irrigated [-]
  - UrbArea:    Catchment fraction classified as urban [-]
  - Qobs:       Mean annual observed streamflow during 1979–2005 [mm/yr]
  - PCPC:       Mean annual CPC precipitation for the period with observed
                streamflow data [mm/yr]
  - PET:        Mean annual potential evaporation during 1979–2005 for the
                period with observed streamflow data [mm/yr]
  - ResInfl:    Ratio of total reservoir capacity to  mean annual observed
                streamflow [-]

The Catchments shapefile also includes the following two fields which are only
provided for the donor catchments:
  - CPC_CAL_01: Daily AOF score using CPC precipitation for the calibration period
  - CPC_VAL_01: Daily AOF score using CPC precipitation for the validation period

The Catchments shapefile also provides, for the evaluation catchments, scores
for 21 performance metrics for two forcing datasets and two parameter sets
(spatially-uniform and regionalized based on the ten most similar donors). Note
that, since the evaluation catchments were not used for deriving the parameters,
we do not distinguish between calibration and validation periods. The fields
were named according to <forcing>_<parameters>_<metric>, where:
        <forcing> denotes the used forcing dataset:
            CPC: CPC precipitation
            WFD: WFDEI precipitation
        <parameters> denotes the parameters:
            UNI: Spatially-uniform parameters
            REG: Regionalized parameters
        <metric> denotes the performance metric:
            01: AOF daily
            02: NSE daily
            03: NSE 5-day
            04: NSE monthly
            05: NSE log-transformed daily
            06: NSE log-transformed 5-day
            07: NSE log-transformed monthly
            08: KGE daily
            09: KGE 5-day
            10: KGE monthly
            11: KGE log-transformed daily
            12: KGE log-transformed 5-day
            13: KGE log-transformed monthly
            14: R² daily
            15: R² 5-day
            16: R² monthly
            17: R² log-transformed daily
            18: R² log-transformed 5-day
            19: R² log-transformed monthly
            20: B
            21: B'
For example, the field WFD_UNI_15 represents the R² score computed from 5-day
flows based on HBV using spatially-uniform parameters and WFDEI precipitation.