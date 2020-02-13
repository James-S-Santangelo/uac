## Metadata

Metadata is only shown for clean/processed datasets that are part of the in
[data-clean](data-clean/), assuming the directory currently exists and the
datasets have been generated. If the datasets have not been generated, run Step
1 in [main.R](scripts/main.R)

### [AllCities_AllPlants.csv](data-clean/AllCities_AllPlants.csv)

| Column | Description | Type |
|--------|-------------|------|
| City  | City from which plants were collected | Factor |
| Population | Population from which plant was collected | Integer |
| Plant | Plant ID | Integer |
| Transect | Transect from which plants were collected. Only relevant for 'Toronto' | Factor |
| Lat.pop | Latitude of population in decimal degrees | Float |
| Long.pop | Longitude of population in decimal degrees | Float |
| HCN_Result | Presence (1) or absence (0) of HCN | Integer |
| Locus.Li | Presence (1) or absence (0) of Linamarase | Integer |
| Locus.Ac | Presence (1) or absence (0) of cyanogenic glucosides | Integer |
| Dmg.Avg | Mean percent leaf area lost (i.e. herbivory), where applicable | Float |
| Distance | Distance (km) of plant to urban centre | Float |


### [AllCities_AllPopulations.csv](data-clean/AllCities_AllPopulations.csv)

| Column | Description | Type |
|--------|-------------|------|
| City  | City from which plants were collected | Factor |
| Transect | Transect from which plants were collected. Only relevant for 'Toronto' | Factor |
| Population Population along urbanization gradient | Integer |
| Distance | Distance (km) of population to urban centre | Float |
| std_distance | Standardized distance of population to city centre | Float |
| n_HCN | Total number of plants in population phenotyped for HCN | Integer |
| sumHCN | Number of cyanogenic plants in population | Integer |
| freqHCN | Frequency of HCN in population. `sumHCN / n_HCN` | Float |
| n_Ac | Total number of plants in population phenotyped for cyanogenic glucosides | Integer |
| sum_Ac | Number of plants in population producing glucosides | Integer |
| freqAc_marker | Frequency of Ac+ in population. `sum_Ac / n_Ac` | Float |
| sum_acac | Total number of Ac– plants in population | Integer |
| freq_acac | Frequency of Ac– plants in population.  `sum_acac / n_Ac` | Float |
| acHWE | Hardy-Weinberg frequency of 'ac' allele. `sqrt(freq_acac)` | Float |
| AcHWE | Hardy-Weinberg frequency of 'Ac' allele. `1 - acHWE` | Float |
| n_Li | Total number of plants in population phenotyped for Linamarase | Integer |
| sum_Li | Number of plants in population producing linamarase | Integer |
| freqLi_marker | Frequency of Li+ in population. `sum_Li / n_Li` | Float |
| sum_lili | Total number of Li– plants in population  | Integer |
| freq_lili | Frequency of Li– plants in population.  `sum_lili / n_Li` | Float |
| liHWE | Hardy-Weinberg frequency of 'li' allele. `sqrt(freq_lili)` | Float |
| LiHWE | Hardy-Weinberg frequency of 'Li' allele. `1 - liHWE` | Float |
| Distance_squared | Squared distance (km) of population to urban centre | Float |
| std_distance_squared | squared standardized distance to urban centre | Float |


### [citySummaryData.csv](data-clean/citySummaryData.csv)

| Column | Description | Type |
|--------|-------------|------|
| City| City from which plants were collected | Factor |
| betaLin_freqHCN | Linear term of best fit model of `freqHCN` against `std_distance` | Float |
| betaQuad_freqHCN | Quadratic term of best fit model (if applicable) of `freqHCN` against `std_distance`  and `std_distance_squared` | Float |
| betaLinOnly | Linear term from model of `freqHCN` against `std_distance`, regardless of best fit | Float |
| betaLog |  Change in log-odds of being cyanogenic from logistic regression of HCN status (0 or 1) against `std_distance` | Float |
| stdErrLog  | Standard error of `betaLog` | Float |
| zLog |  z-statistic associated with `betaLog` | Float |
| pvalLog | P-value associated with `betaLog` | Float |
| freqHCN | Frequency of HCN in population. `sumHCN / n_HCN` | Float |
| AcHWE | Hardy-Weinberg frequency of 'Ac' allele. `1 - acHWE` | Float |
| LiHWE | Hardy-Weinberg frequency of 'Li' allele. `1 - liHWE` | Float |
| Latitude | Latitude of city in decimal degrees | Float |
| Longitude | Longitude of city in decimal degrees | Float |
| annualAI | Annual aridity index of city | Float |
| monthlyPET | Mean monthly potential evapotranspiration of city | Float |
| annualPET | Annual potential evapotranspiration of city  | Float |
| monthlyPrecip | Mean monthly summer precipitation of city  | Float |
| mwtBio | Mean minimum winter temperature (°C) of city | Float |
| mstBio | Mean maximum summer temperature (°C) of city | Float |
| smd | Soil moisture deficit | Float |
| snow_depth | Mean winter snow depth (cm) | Float |
| snowfall | Mean winter snowfall (cm) | Float |
| mwtWea | Mean minimum winter temperature (°C) from weather station data | Float |
| daysNegNoSnow | Mean number of days < 0°C with no snow cover | Float |


### [haplotypeData.csv](data-clean/haplotypeData.csv)

| Column | Description | Type |
|--------|-------------|------|
| City  | City from which plants were collected | Factor |
| Population | Population from which plant was collected | Integer |
| Plant | Plant ID | Integer |
| HCN_Result | Presence (1) or absence (0) of HCN | Integer |
| ac1.3 | Presence (1) or absence (0) of PCR product for 1.3 kb deletion at *Ac* locus | Integer |
| ac2.2 | Presence (1) or absence (0) of PCR product for 2.2 kb deletion at *Ac* locus | Integer |
| ac5.3 | Presence (1) or absence (0) of PCR product for 5.3 kb deletion at *Ac* locus | Integer |
| li4.0 | Presence (1) or absence (0) of PCR product for 4.0 kb deletion at *Li* locus | Integer |
| ac7.0 | Presence (1) or absence (0) of PCR product for 7.0 kb deletion at *Li* locus | Integer |
| ac9.3 | Presence (1) or absence (0) of PCR product for 9.3 kb deletion at *Li* locus | Integer |
| haplotype_Ac | Concatenated haplotype of plant at the *Ac* Locus | Factor |
| haplotype_Li | Concatenated haplotype of plant at the *Li* Locus | Factor |
| haplotype_all | Concatenated haplotype of plant at both loci | Factor |
| Distance | Distance (km) of population to urban centre | Float |
| Habitat | One of 'urban' or 'rural' | Factor







