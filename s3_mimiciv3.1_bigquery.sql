WITH sep AS (
    SELECT *,
    LEAST(se.suspected_infection_time, se.sofa_time) AS sepsis_onset
    FROM `physionet-data.mimiciv_3_1_derived.sepsis3` AS se
),
s3 AS (
    SELECT 
        icu.*,
        sep.sepsis_onset,
        sep.sofa_score,
        sep.sepsis3,
        CASE
            WHEN icu.hospital_expire_flag = 1 THEN 1
            WHEN icu.hospital_expire_flag = 0 AND icu.dod <= icu.dischtime THEN 1
            WHEN icu.hospital_expire_flag = 0 AND icu.dod > icu.dischtime THEN 0
            WHEN icu.dod IS NULL THEN 0
            ELSE NULL
        END AS hospital_outcome,
        CASE 
            WHEN icu.dod <= icu.icu_outtime THEN 1
            ELSE 0 
        END AS icu_outcome
    FROM `physionet-data.mimiciv_3_1_derived.icustay_detail` AS icu
    INNER JOIN sep 
        ON icu.stay_id = sep.stay_id
    WHERE sep.sepsis3 = TRUE 
        AND icu.first_icu_stay = TRUE
        AND sep.sepsis_onset >= DATETIME_SUB(icu.icu_intime, INTERVAL 6 HOUR)
        AND sep.sepsis_onset <= DATETIME_ADD(icu.icu_intime, INTERVAL 1 DAY)
        AND icu.admission_age >= 18
        AND icu.los_icu >= 1
),
s3_1 AS (
SELECT 
    *,
    CASE 
        WHEN hospital_outcome = 1 THEN DATETIME_DIFF(dischtime, icu_intime, DAY)
        WHEN hospital_outcome = 0 THEN 360
        ELSE NULL 
    END AS survival_time_crude
FROM s3
),
s3_2 AS (
    SELECT *,
    CASE 
        WHEN survival_time_crude >= 28 THEN 28
        ELSE survival_time_crude 
    END AS survival_time_28d,
    CASE 
        WHEN survival_time_crude >= 7 THEN 7
        ELSE survival_time_crude 
    END AS survival_time_7d,
    CASE 
        WHEN survival_time_crude >= 90 THEN 90
        ELSE survival_time_crude 
    END AS survival_time_90d
    FROM s3_1
),
lab_cbc AS (
    SELECT
        ie.stay_id,
        AVG(hematocrit) AS hematocrit_mean,
        AVG(hemoglobin) AS hemoglobin_mean,
        AVG(platelet) AS platelets_mean,
        AVG(wbc) AS wbc_mean
    FROM
        s3_2 AS ie
        LEFT JOIN `physionet-data.mimiciv_3_1_derived.complete_blood_count` AS le 
            ON le.subject_id = ie.subject_id
            AND le.charttime >= DATETIME_SUB(ie.icu_intime, INTERVAL 6 HOUR)
            AND le.charttime <= DATETIME_ADD(ie.icu_intime, INTERVAL 1 DAY)
    GROUP BY
        ie.stay_id
),
lab_chem AS (
    SELECT
        ie.stay_id,
        AVG(albumin) AS albumin_mean,
        AVG(globulin) AS globulin_mean,
        AVG(total_protein) AS total_protein_mean,
        AVG(aniongap) AS aniongap_mean,
        AVG(bicarbonate) AS bicarbonate_mean,
        AVG(bun) AS bun_mean,
        AVG(calcium) AS calcium_mean,
        AVG(chloride) AS chloride_mean,
        AVG(creatinine) AS creatinine_mean,
        AVG(glucose) AS glucose_mean,
        AVG(sodium) AS sodium_mean,
        AVG(potassium) AS potassium_mean
    FROM
        s3_2 AS ie
        LEFT JOIN `physionet-data.mimiciv_3_1_derived.chemistry` AS le 
            ON le.subject_id = ie.subject_id
            AND le.charttime >= DATETIME_SUB(ie.icu_intime, INTERVAL 6 HOUR)
            AND le.charttime <= DATETIME_ADD(ie.icu_intime, INTERVAL 1 DAY)
    GROUP BY
        ie.stay_id
),
lab_diff AS (
    SELECT
        ie.stay_id,
        AVG(basophils_abs) AS abs_basophils_mean,
        AVG(eosinophils_abs) AS abs_eosinophils_mean,
        AVG(lymphocytes_abs) AS abs_lymphocytes_mean,
        AVG(monocytes_abs) AS abs_monocytes_mean,
        AVG(neutrophils_abs) AS abs_neutrophils_mean,
        AVG(atypical_lymphocytes) AS atyps_mean,
        AVG(bands) AS bands_mean,
        AVG(immature_granulocytes) AS imm_granulocytes_mean,
        AVG(metamyelocytes) AS metas_mean,
        AVG(nrbc) AS nrbc_mean
    FROM
        s3_2 AS ie
        LEFT JOIN `physionet-data.mimiciv_3_1_derived.blood_differential` AS le 
            ON le.subject_id = ie.subject_id
            AND le.charttime >= DATETIME_SUB(ie.icu_intime, INTERVAL 6 HOUR)
            AND le.charttime <= DATETIME_ADD(ie.icu_intime, INTERVAL 1 DAY)
    GROUP BY
        ie.stay_id
),
lab_coag AS (
    SELECT
        ie.stay_id,
        AVG(d_dimer) AS d_dimer_mean,
        AVG(fibrinogen) AS fibrinogen_mean,
        AVG(thrombin) AS thrombin_mean,
        AVG(inr) AS inr_mean,
        AVG(pt) AS pt_mean,
        AVG(ptt) AS ptt_mean
    FROM
        s3_2 AS ie
        LEFT JOIN `physionet-data.mimiciv_3_1_derived.coagulation` AS le 
            ON le.subject_id = ie.subject_id
            AND le.charttime >= DATETIME_SUB(ie.icu_intime, INTERVAL 6 HOUR)
            AND le.charttime <= DATETIME_ADD(ie.icu_intime, INTERVAL 1 DAY)
    GROUP BY
        ie.stay_id
),
lab_enz AS (
    SELECT
        ie.stay_id,
        AVG(alt) AS alt_mean,
        AVG(alp) AS alp_mean,
        AVG(ast) AS ast_mean,
        AVG(amylase) AS amylase_mean,
        AVG(bilirubin_total) AS bilirubin_total_mean,
        AVG(bilirubin_direct) AS bilirubin_direct_mean,
        AVG(bilirubin_indirect) AS bilirubin_indirect_mean,
        AVG(ck_cpk) AS ck_cpk_mean,
        AVG(ck_mb) AS ck_mb_mean,
        AVG(ggt) AS ggt_mean,
        AVG(ld_ldh) AS ld_ldh_mean
    FROM
        s3_2 AS ie
        LEFT JOIN `physionet-data.mimiciv_3_1_derived.enzyme` AS le 
            ON le.subject_id = ie.subject_id
            AND le.charttime >= DATETIME_SUB(ie.icu_intime, INTERVAL 6 HOUR)
            AND le.charttime <= DATETIME_ADD(ie.icu_intime, INTERVAL 1 DAY)
    GROUP BY
        ie.stay_id
),

lab_bg AS (
    SELECT
        ie.stay_id,
        AVG(bg.lactate) AS lactate_mean,
        AVG(bg.ph) AS ph_mean,
        AVG(bg.so2) AS so2_mean,
        AVG(bg.po2) AS po2_mean,
        AVG(bg.pco2) AS pco2_mean,
        AVG(bg.aado2) AS aado2_mean,
        AVG(bg.aado2_calc) AS aado2_calc_mean,
        AVG(bg.pao2fio2ratio) AS pao2fio2ratio_mean,
        AVG(bg.baseexcess) AS baseexcess_mean,
        AVG(bg.totalco2) AS totalco2_mean
    FROM
        s3_2 AS ie
        LEFT JOIN `physionet-data.mimiciv_3_1_derived.bg` AS bg 
            ON bg.subject_id = ie.subject_id
            AND bg.charttime >= DATETIME_SUB(ie.icu_intime, INTERVAL 6 HOUR)
            AND bg.charttime <= DATETIME_ADD(ie.icu_intime, INTERVAL 1 DAY)
    GROUP BY
        ie.stay_id
),
s3_3 AS (
    SELECT
        s3_2.*, 
        lab_cbc.hematocrit_mean, lab_cbc.hemoglobin_mean, lab_cbc.platelets_mean, lab_cbc.wbc_mean,
        lab_chem.albumin_mean, lab_chem.globulin_mean, lab_chem.total_protein_mean, lab_chem.aniongap_mean,
        lab_chem.bicarbonate_mean, lab_chem.bun_mean, lab_chem.calcium_mean, lab_chem.chloride_mean,
        lab_chem.creatinine_mean, lab_chem.glucose_mean, lab_chem.sodium_mean, lab_chem.potassium_mean,
        lab_diff.abs_basophils_mean, lab_diff.abs_eosinophils_mean, lab_diff.abs_lymphocytes_mean,
        lab_diff.abs_monocytes_mean, lab_diff.abs_neutrophils_mean, lab_diff.atyps_mean, lab_diff.bands_mean,
        lab_diff.imm_granulocytes_mean, lab_diff.metas_mean, lab_diff.nrbc_mean, lab_coag.d_dimer_mean,
        lab_coag.fibrinogen_mean, lab_coag.thrombin_mean, lab_coag.inr_mean, lab_coag.pt_mean, lab_coag.ptt_mean,
        lab_enz.alt_mean, lab_enz.alp_mean, lab_enz.ast_mean, lab_enz.amylase_mean, lab_enz.bilirubin_total_mean,
        lab_enz.bilirubin_direct_mean, lab_enz.bilirubin_indirect_mean, lab_enz.ck_cpk_mean, lab_enz.ck_mb_mean,
        lab_enz.ggt_mean, lab_enz.ld_ldh_mean
    FROM s3_2
    LEFT JOIN lab_cbc ON s3_2.stay_id = lab_cbc.stay_id
    LEFT JOIN lab_chem ON s3_2.stay_id = lab_chem.stay_id
    LEFT JOIN lab_diff ON s3_2.stay_id = lab_diff.stay_id
    LEFT JOIN lab_coag ON s3_2.stay_id = lab_coag.stay_id
    LEFT JOIN lab_enz ON s3_2.stay_id = lab_enz.stay_id
),
s3_3_1 AS (
    SELECT 
        s3_3.*,
        lab_bg.lactate_mean, lab_bg.ph_mean, lab_bg.so2_mean, lab_bg.po2_mean, lab_bg.pco2_mean, lab_bg.aado2_mean,
        lab_bg.aado2_calc_mean, lab_bg.pao2fio2ratio_mean, lab_bg.baseexcess_mean, lab_bg.totalco2_mean
    FROM s3_3
    LEFT JOIN lab_bg ON s3_3.stay_id = lab_bg.stay_id
),
s3_4 AS (
SELECT 
    s3_3_1.*,
    gcs.gcs_min
FROM s3_3_1
LEFT JOIN `physionet-data.mimiciv_3_1_derived.first_day_gcs` AS gcs
    ON s3_3_1.stay_id = gcs.stay_id
),
s3_5 AS(
SELECT
    s3_4.*,
    vital.heart_rate_mean,
    vital.sbp_mean,
    vital.dbp_mean,
    vital.mbp_mean,
    vital.resp_rate_mean,
    vital.temperature_mean,
    vital.spo2_mean
FROM s3_4
LEFT JOIN `physionet-data.mimiciv_3_1_derived.first_day_vitalsign` AS vital
    ON s3_4.stay_id = vital.stay_id
),
s3_6 AS (
SELECT 
    s3_5.*,
    hei.height
FROM s3_5
LEFT JOIN `physionet-data.mimiciv_3_1_derived.first_day_height` AS hei
    ON s3_5.stay_id = hei.stay_id
),
s3_7 AS (
SELECT 
    s3_6.*,
    wei.weight
FROM s3_6
LEFT JOIN `physionet-data.mimiciv_3_1_derived.first_day_weight` AS wei
    ON s3_6.stay_id = wei.stay_id
),
s3_8 AS (
SELECT 
    s3_7.*,
    urin.urineoutput
FROM s3_7
LEFT JOIN `physionet-data.mimiciv_3_1_derived.first_day_urine_output` AS urin
    ON s3_7.stay_id = urin.stay_id
),
s3_9 AS (
SELECT 
    s3_8.*,
    rrt.dialysis_active,
    rrt.dialysis_type
FROM s3_8
LEFT JOIN `physionet-data.mimiciv_3_1_derived.first_day_rrt` AS rrt
    ON s3_8.stay_id = rrt.stay_id
),
s3_10 AS(
SELECT 
    s3_9.*,
    aps.apsiii 
FROM s3_9
LEFT JOIN `physionet-data.mimiciv_3_1_derived.apsiii` AS aps 
    ON s3_9.stay_id = aps.stay_id
),
s3_11 AS(
SELECT 
    s3_10.*,
    lod.lods
FROM s3_10
LEFT JOIN `physionet-data.mimiciv_3_1_derived.lods` AS lod
    ON s3_10.stay_id = lod.stay_id
),
s3_12 AS(
SELECT 
    s3_11.*,
    oa.oasis
FROM s3_11
LEFT JOIN `physionet-data.mimiciv_3_1_derived.oasis` AS oa
    ON s3_11.stay_id = oa.stay_id
),
-- SIRS (Systemic Inflammatory Response Syndrome) score calculation
-- This query extracts the SIRS criteria
-- The score is calculated on the first day of each ICU patients' stay
sirs_scorecomp AS (
  SELECT
    ie.stay_id,
    v.temperature_min,
    v.temperature_max,
    v.heart_rate_max,
    v.resp_rate_max,
    bg.pco2_min AS paco2_min,
    l.wbc_min,
    l.wbc_max,
    l.bands_max
  FROM `physionet-data.mimiciv_3_1_icu.icustays` AS ie
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.first_day_bg_art` AS bg
    ON ie.stay_id = bg.stay_id
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.first_day_vitalsign` AS v
    ON ie.stay_id = v.stay_id
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.first_day_lab` AS l
    ON ie.stay_id = l.stay_id
),
sirs_scorecalc AS (
  -- Calculate the final score
  -- note that if the underlying data is missing, the component is null
  -- eventually these are treated as 0 (normal), but knowing when
  -- data is missing is useful for debugging
  SELECT
    stay_id,
    CASE
      WHEN temperature_min < 36.0
      THEN 1
      WHEN temperature_max > 38.0
      THEN 1
      WHEN temperature_min IS NULL
      THEN NULL
      ELSE 0
    END AS temp_score,
    CASE
      WHEN heart_rate_max > 90.0
      THEN 1
      WHEN heart_rate_max IS NULL
      THEN NULL
      ELSE 0
    END AS heart_rate_score,
    CASE
      WHEN resp_rate_max > 20.0
      THEN 1
      WHEN paco2_min < 32.0
      THEN 1
      WHEN COALESCE(resp_rate_max, paco2_min) IS NULL
      THEN NULL
      ELSE 0
    END AS resp_score,
    CASE
      WHEN wbc_min < 4.0
      THEN 1
      WHEN wbc_max > 12.0
      THEN 1
      WHEN bands_max > 10
      THEN 1
      -- > 10% immature neutrophils (band forms)
      WHEN COALESCE(wbc_min, bands_max) IS NULL
      THEN NULL
      ELSE 0
    END AS wbc_score
  FROM sirs_scorecomp
),
sirs_table AS (
  SELECT
    ie.stay_id,
    -- Combine all the scores to get SIRS
    -- Impute 0 if the score is missing
    COALESCE(temp_score, 0) + COALESCE(heart_rate_score, 0) + COALESCE(resp_score, 0) + COALESCE(wbc_score, 0) AS sirs,
    temp_score,
    heart_rate_score,
    resp_score,
    wbc_score
  FROM `physionet-data.mimiciv_3_1_icu.icustays` AS ie
  LEFT JOIN sirs_scorecalc AS s
    ON ie.stay_id = s.stay_id
),
s3_14 AS(
SELECT 
    s3_12.*,
    si.sirs
FROM s3_12
LEFT JOIN sirs_table AS si
    ON s3_12.stay_id = si.stay_id
),
s3_15 AS(
SELECT 
    s3_14.*,
    cha.myocardial_infarct,
    cha.congestive_heart_failure,
    cha.peripheral_vascular_disease,
    cha.cerebrovascular_disease,
    cha.dementia,
    cha.chronic_pulmonary_disease,
    cha.rheumatic_disease,
    cha.peptic_ulcer_disease,
    cha.mild_liver_disease,
    cha.diabetes_without_cc,
    cha.diabetes_with_cc,
    cha.paraplegia,
    cha.renal_disease,
    cha.malignant_cancer,
    cha.severe_liver_disease,
    cha.metastatic_solid_tumor,
    cha.aids,
    cha.charlson_comorbidity_index
FROM s3_14
LEFT JOIN `physionet-data.mimiciv_3_1_derived.charlson` AS cha
    ON s3_14.hadm_id = cha.hadm_id
),
diag AS (
  SELECT
    hadm_id,
    CASE WHEN icd_version = 9 THEN icd_code ELSE NULL END AS icd9_code,
    CASE WHEN icd_version = 10 THEN icd_code ELSE NULL END AS icd10_code
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
),
exclu AS (
  SELECT
    hadm_id,
    MAX(
      CASE
        WHEN icd10_code IN ('D51', 'D52', 'D59') OR
          SUBSTR(icd10_code, 1, 3) BETWEEN 'D60' AND 'D64' OR
          SUBSTR(icd10_code, 1, 3) BETWEEN 'D70' AND 'D89' OR
          SUBSTR(icd9_code, 1, 4) BETWEEN '2810' AND '2812' OR
          icd9_code = '2830' OR
          icd9_code LIKE '284%' OR
          icd9_code = '2851' OR
          icd9_code IN ('28521', '28522', '28529', '2853') OR
          icd9_code = '2882' OR
          icd9_code = '2884' OR
          icd9_code = '2894' OR
          icd9_code IN ('28950', '28951', '28959')
        THEN 1
        ELSE 0
      END
    ) AS hematology_disease,
    MAX(
      CASE
        WHEN icd10_code IN ('Z510', 'Z511', 'Z5111', 'Z5112') OR
          icd9_code IN ('V580', 'V581', 'V5811', 'V5812')
        THEN 1
        ELSE 0
      END
    ) AS anti_cancer_therapy,
    MAX(
      CASE
        WHEN icd10_code IN ('Z795', 'Z7951', 'Z7952') OR
          icd9_code = 'V5865'
        THEN 1
        ELSE 0
      END
    ) AS longterm_steroid,
    MAX(
      CASE
        WHEN SUBSTR(icd10_code, 1, 3) = 'Z94' OR
         SUBSTR(icd9_code, 1, 3) = 'V42'
        THEN 1
        ELSE 0
      END
    ) AS transplant
  FROM diag
  GROUP BY hadm_id
),
s3_16 AS (
    SELECT 
        s3_15.*,
        COALESCE(exclu.hematology_disease, 0) AS hematology_disease,
        COALESCE(exclu.anti_cancer_therapy, 0) AS anti_cancer_therapy,
        COALESCE(exclu.longterm_steroid, 0) AS longterm_steroid,
        COALESCE(exclu.transplant, 0) AS transplant
    FROM s3_15
    LEFT JOIN exclu 
        ON s3_15.hadm_id = exclu.hadm_id
),
SELECT * FROM s3_16;


