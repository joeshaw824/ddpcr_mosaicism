# ddpcr_mosaicism

Analysis of ddPCR data for confirmation of mosaic DNA variants detected by high-depth NGS sequencing.

Taqman assays (ThermoFisher) are designed for specific mosaic variants and optimised by performing ddPCR on patient, normal and NTC (no template control) samples at a range of annealing temperatures (65-55 degrees Celsius).

Assay ID (assay_id) is used as a consistent identifier of a specific assay design. The assay ID is a unique 7 digit alphanumeric string which is specific to a Taqman assay design and is assigned by ThermoFisher when the assay is ordered. The assay name (assay_name) is the name given to an assay design by the user, and may not be unique.
