# Fermi-LAT Catalog Resources

This directory contains local copies of the Fermi-LAT catalog FITS files used by the Fermi catalog workflow.

## Files

| Local file             | Catalog  | Description                                                   | Original FSSC file      |
| ---------------------- | -------- | ------------------------------------------------------------- | ----------------------- |
| `4FGL_DR4.fit`         | 4FGL-DR4 | Fermi-LAT 14-year Source Catalog                              | `gll_psc_v35.fit`       |
| `4LAC_DR3_H.fits`      | 4LAC-DR3 | Fourth LAT AGN Catalog, Data Release 3, high-latitude sources | `table-4LAC-DR3-h.fits` |
| `4LAC_DR3_L.fits`      | 4LAC-DR3 | Fourth LAT AGN Catalog, Data Release 3, low-latitude sources  | `table-4LAC-DR3-l.fits` |
| `4LAC_DR3_merged.fits` | 4LAC-DR3 | Merged high- and low-latitude 4LAC-DR3 catalog with an additional `LAC_sample` column indicating the source catalog | derived locally |
| `3FHL.fit`             | 3FHL     | Third Fermi-LAT Catalog of High-Energy Sources                | `gll_psch_v13.fit`      |

## Source

The files were obtained from the official Fermi Science Support Center (FSSC) catalog data pages:

* 4FGL-DR4: https://fermi.gsfc.nasa.gov/ssc/data/access/lat/14yr_catalog/
* 4LAC-DR3: https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4LACDR3/
* 3FHL: https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3FHL/

## Notes

The local filenames are normalized for use inside this repository. The corresponding original FSSC filenames are listed above to keep the catalog provenance explicit.

`4LAC_DR3_merged.fits` is generated locally by vertically merging the high- and low-latitude 4LAC-DR3 FITS tables. An additional `LAC_sample` column records whether each source originates from the `high_lat` or `low_lat` sample. The merged file is not an independently distributed FSSC catalog product.
