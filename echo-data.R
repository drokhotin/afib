library(tidyr)
library(dplyr)
library(jsonlite)
library(ggplot2)

echodf <- fromJSON(txt = "/home/okhotin/Dropbox/R/masters/echo.json")

echo <- echodf[,c("rec_id", "pat_id")]
colnames(echo) <- c("record_id", "patient_id")

echo <- echo %>% 
  mutate(id = row_number()) ## add unique id

echo$weight <- as.numeric(echodf$value$studydata$HOWeight)
echo$height <- as.numeric(echodf$value$studydata$HOHeight)
echo$bsa <- as.numeric(echodf$value$studydata$bsa)
echo$date <- as.Date(echodf$value$studydata$ex_date, format="%Y%m%d")

# Averaging BSA, weight and height, based on assumption that these
# parameters are relatively stable across the follow-up (checked for
# patients from the study -- difference of BSA no more then 0.25 from averaged)

echo <- echo %>%
  group_by(patient_id) %>%
  mutate(bsa = round(mean(bsa, na.rm=T),2),
         weight = round(mean(weight, na.rm=T)),
         height = round(mean(height, na.rm=T))) %>%
  ungroup()

# EDV index

# Teicholz EDV estimate is biased so it was corrected
# based on linear regression derived from studies where 
# both measurements values were present
  
p_edv_bi <- echodf$value$measurements$`Simpson BP`$`EDV`$value
p_edv_4 <- echodf$value$measurements$`Simpson BP`$`EDV4`$value
p_edv_s <- echodf$value$measurements$`Simpson SP`$`EDV`$value
p_edv_t <- echodf$value$measurements$`Teichholz(2D)`$`EDV`$value

edv <- data.frame(edv_bi=as.numeric(p_edv_bi), edv_4=as.numeric(p_edv_4),
                  edv_s=as.numeric(p_edv_s), edv_t=as.numeric(p_edv_t))

slope <- summary(lm(data= edv, edv_bi ~ edv_t))$coef[2][1]
intercept <- summary(lm(data= edv, edv_bi ~ edv_t))$coef[1][1]

edv$edv_bit <- intercept + slope * edv$edv_t

echo$edv <- ifelse(!is.na(edv$edv_bi), edv$edv_bi,
                ifelse(!is.na(edv$edv_4), edv$edv_4, 
                       ifelse(!is.na(edv$edv_s), edv$edv_s, edv$edv_bit)
                ))

echo$edvi <- echo$edv/echo$bsa # calculating EDV index based on BSA



## Ejection fraction

p_ef_bi <- echodf$value$measurements$`Simpson BP`$EF$value
p_ef_4 <- echodf$value$measurements$`Simpson BP`$EF4$value
p_ef_s <- echodf$value$measurements$`Simpson SP`$EF$value

p_ef <- ifelse(!is.na(p_ef_bi), p_ef_bi,
                ifelse(!is.na(p_ef_4), p_ef_4, p_ef_s))

echo$ef <- as.numeric(p_ef)


## LV Mass

p_lvm_al <- as.numeric(echodf$value$measurements$`LV Mass A-L`$`LV Mass`$value)
p_lvm_te <- as.numeric(echodf$value$measurements$`Teichholz(2D)`$`LV Mass`$value)

lvm <- data.frame(al=p_lvm_al, te=p_lvm_te)

slope <- summary(lm(data= lvm, al ~ te))$coef[2][1]
intercept <- summary(lm(data= lvm, al ~ te))$coef[1][1]

lvm$te_c <- intercept + slope * lvm$te

lvm$lvm <- ifelse(!is.na(lvm$al), lvm$al, lvm$te_c)

echo$lvmi <- lvm$lvm / echo$bsa

##


p_lav_s <- echodf$value$measurements$`LA Vol (Simp)`$`LA Vol (Simp)`$value
p_lav_s4 <- echodf$value$measurements$`LA Vol (Simp)`$`4CH`$value
p_lav_s2 <- echodf$value$measurements$`LA Vol (Simp)`$`2CH`$value
p_lav_AL <- echodf$value$measurements$`LA Vol (A-L)`$`LA Vol (A-L)`$value
p_lav_LAD <- echodf$value$measurements$`AV/LA`$`LAs diam`$value


lav <- data.frame(as.numeric(p_lav_s), as.numeric(p_lav_s4), 
                   as.numeric(p_lav_s2), as.numeric(p_lav_AL),
                   as.numeric(p_lav_LAD))
colnames(lav) <- c('s', 's4', 's2', 'al', 'lad')

lav$lav <- ifelse(!is.na(lav$s), lav$s, 
                    ifelse(!is.na(lav$s4), lav$s4, 
                           ifelse(!is.na(lav$s2), lav$s2, lav$al)
                    ))

echo$lavi <- lav$lav / echo$bsa

### Removing results of miscalculating indices (by dividing them by NA BSA)

echo$lavi[is.nan(echo$lavi)] <- NA
echo$lvmi[is.nan(echo$lvmi)] <- NA
echo$weight[is.nan(echo$weight)] <- NA
echo$height[is.nan(echo$height)] <- NA
echo$bsa[is.nan(echo$bsa)] <- NA


######### Single parameter for every patients
## From the first study where this parameter is encountered



echo_single <- echo %>%
  group_by(patient_id) %>%
  mutate(
    study_n = row_number(),
    edvi=first(na.omit(edvi)),
    ef=first(na.omit(ef)),
    lvmi=first(na.omit(lvmi)),
    lavi=first(na.omit(lavi))
  ) %>%
  filter(study_n==1)  %>%
  select(patient_id, weight, height, bsa,  edvi, ef, lvmi, lavi)




