setGeneric("DispX_new", function(ODD,Omega,center, BD_params, LL,Method, epsilon=c(0.15,0.03,0.1))
  standardGeneric("DispX_new") )
# Code that calculates/predicts the total human displacement 
setMethod("DispX_new", "ODD", function(ODD,Omega,center, BD_params, LL=F,
                                   Method=list(Np=1,cores=1,cap=-300), epsilon=c(0.15,0.03,0.1)
){
  # Extract 0D parameters & speed up loop
  Params<-FormParams(ODD,list(Np=Method$Np,center=center))
  # Income distribution percentiles & extract income percentile  
  SincN<-seq(10,90,by = 10); Sinc<-ExtractCIndy(ODD,var = paste0("p",SincN,"p100"))
  # Speed-up calculation (through accurate cpu-work distribution) to only values that are not NA
  if(!LL) {notnans<-which(!(is.na(ODD@data[["Population"]]) | is.na(ODD@data[["ISO3C"]]) | is.na(ODD@data[["GDP"]])))
  } else notnans<-which(!(is.na(ODD@data[["Population"]]) | is.na(ODD@data[["ISO3C"]]) | is.na(ODD@data[["GDP"]]) |
                            !ODD@data[["ISO3C"]]%in%ODD@gmax$iso3))
  
  BD_data_present <- ifelse(all(is.na(ODD@data$nBuildings)) , F, T)
  # Calculate non-local linear predictor values
  LP<-GetLP(ODD,Omega,Params,Sinc,notnans)
  # Speed things up a little
  hrange<-grep("hazMean",names(ODD),value = T)
  # Function to predict displacement per gridpoint
  CalcDam<-function(ij){
    iso3c<-ODD@data$ISO3C[ij]
    # Calculate local linear predictor (NOTE: is a vector due to income distribution)
    locallinp<-LP$dGDP$linp[LP$dGDP$ind==LP$iGDP[ij]]*LP$Plinp[ij]*LP$linp[[iso3c]] #LOOSEEND
    #locallinp<-rep(1,10) #reducing parameter space while I'm figuring out the MCMC
    
    # Sample population per income distribution (Assumes 9 percentiles):
    lPopS <- SplitSamplePop(Pop=ODD@data$Population[ij],Method$Np) 
    tPop <-array(0,c(3, Method$Np)) #row 1 = tDisp, #row 2 = tMort, #row 3 = tRem
    tPop[3,]=colSums(lPopS)
    for(h in hrange){
      # for(h in c(1)){
      if(is.na(ODD@data[ij,h])) next
      # Resample population based on who is remaining
      ind<-tPop[3,]>0
      if(h!=hrange[1]) {
        if(sum(ind)==0) break #if no remaining population, skip modelling
        if(length(lPopS[,ind])==0) break #if no remaining population, skip modelling
        #if(sum(ind)>1) sumz<-colSums(lPopS[,ind])
        #else sumz<-sum(lPopS[,ind])
        #lPopS[,!ind]<-0
        lPopS[,ind]<-SplitSamplePop(Pop=tPop[3,ind])
      }
      # Sample hazard Intensity 
      # the uncertainty is too high... so I scale it to get some interpretable results (I know, I'm not really a statistician, I don't even have a degree, I was actually just having a look around the department when they confused me for the interviewee. I didn't have the heart to say anything. You don't hate me as much as I do)
      # I_ij<-rnorm(n = Method$Np,
      #             mean = ODD@data[ij,paste0("hazMean",h)],
      #             sd = ODD@data[ij,paste0("hazSD",h)]/10)
      I_ij<-ODD@data[ij,h]
      
      # Separate into income distributions (as each have 10% of population, order doesn't matter)
      for (s in 1:length(SincN)){
        if(all(lPopS[s,]==0)) next
        # Predict damage at coordinate {i,j} (vector with MC particles)
        Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=Params$Np),Omega)*locallinp[s], error=function(e) NA)
        if(any(is.na(Damage))) print(ij)
        
        D_MortDisp <- D_MortDisp_calc(Damage, Omega) #First row of D_MortDisp is D_Mort, second row is D_Disp
        
        # Accumulate the number of people displaced/deceased, but don't accumulate the remaining population
        tPop[3,ind]<-0
        
        tPop[,ind]<-tPop[,ind] + Fbdam(lPopS[s,ind],D_MortDisp[2,ind], D_MortDisp[1,ind], (1-D_MortDisp[1,ind]-D_MortDisp[2,ind]))
      } 
    }
    #ensure the total displaced, deceased or remaining does not exceed total population
    tPop[tPop>ODD@data$Population[ij]]<-floor(ODD@data$Population[ij])
    
    #if no building destruction data:
    if(!BD_data_present) return(rbind(tPop[1:2,, drop=FALSE], rep(NA, Method$Np))) #return total displacement and mortality, set number of buildings destroyed to NA
    
    #otherwise, sample from the model for the number of buildings destroyed:
    #we take locallinp[5] which corresponds to locallinp for the median GDP
    nBuildings = rep(ODD@data$nBuildings[ij], Method$Np)
    nBD = rep(0, Method$Np)
    for (h in hrange){
      if(is.na(ODD@data[ij,h])) next
      if(h!=hrange[1]) {
        if(all(nBuildings==0)) break #if no remaining buildings, skip modelling LOOSEEND
      }
      
      I_ij<-ODD@data[ij,h]
      Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=Params$Np),Omega)*locallinp[5], error=function(e) NA) #calculate unscaled damage (excluding GDP)
      
      D_BD = plnorm(Damage, meanlog = Omega$Lambda3$nu, sdlog = Omega$Lambda3$omega)
      
      moreBD = fBD(nBuildings, D_BD)
      nBD = nBD + moreBD
      
      nBuildings = nBuildings - moreBD
    }
    return(rbind(tPop[1:2,,drop=FALSE], nBD))
  }
  
  Dam<-array(0,c(nrow(ODD),Method$Np,3)) # Dam[,,1] = Displacement, Dam[,,2] = Mortality, Dam[,,3] = Buildings Destroyed
  
  #Method$cores is equal to AlgoParams$NestedCores (changed in Model file)
  if(Method$cores>1) { Dam[notnans,,]<-aperm(simplify2array(mclapply(X = as.data.frame(notnans), FUN = CalcDam,mc.cores = Method$cores)), perm=c(3,2,1))
  } else  Dam[notnans,,]<- aperm(simplify2array(lapply(X = as.data.frame(notnans), FUN = CalcDam)), perm=c(3,2,1))
  
  # return(Disp)
  
  # If the IDMC estimate is foreseen to be a lower or upper bound, or a generally poor estimate
  # for(c in ODD@gmax$iso3){
  #   ind<-ODD@data$ISO3C==c  & !is.na(ODD@data$ISO3C)
  #   Disp[ind,]%<>%qualifierDisp(qualifier = ODD@gmax$qualifier[ODD@gmax$iso3==c],mu = Omega$mu)
  # }
  
  funcy<-function(i,LLout=T, epsilon=AlgoParams$epsilon_min, kernel=AlgoParams$kernel, cap=AlgoParams$cap) {
    tmp<-data.frame(iso3=ODD@data$ISO3C,IDPs=Dam[,i,1], mort=Dam[,i,2], nBD=Dam[,i,3]) %>% 
      group_by(iso3) %>% summarise(disp_predictor=floor(sum(IDPs,na.rm = T)), 
                                   mort_predictor=floor(sum(mort,na.rm = T)),
                                   nBD_predictor=floor(sum(nBD,na.rm = T)),
                                   .groups = 'drop_last')
    tmp<-tmp[!is.na(tmp$iso3) & tmp$iso3%in%ODD@gmax$iso3,]
    tmp%<>%merge(ODD@gmax,by="iso3")%>%arrange(desc(gmax))
    #print(paste(tmp$nBD_predictor, tmp$buildDestroyed))
    #print(tmp)
    if(LLout) {
      return(LL_IDP(tmp, epsilon,  kernel, cap))
    }
    return(tmp)
  }
  if (LL == F & Method$Np == 1){
    ODD@predictDisp<-funcy(1,LLout=F) 
    
    return(ODD)
  }
  
  outer<-vapply(1:Method$Np,funcy,numeric(length(unique(ODD@gmax$iso3))), epsilon=epsilon)
  outer[outer < (-745)]<- -745
  
  # Find thebest fit solution
  if(length(unique(ODD@gmax$iso3))>1) {
    if(LL)  return(log(rowMeans(exp(outer),na.rm=T)))
    MLE<-which.max(log(colSums(exp(outer),na.rm=T)))
  }  else {
    #return(outer) #SMC-CHANGE
    if(LL)  return(log(mean(exp(outer),na.rm=T)))
    MLE<-which.max(log(exp(outer)))
  }
  if(Method$Np == 1){
    MLE=1
  }
  # Save into ODD object
  # ODD@data$Disp<-Disp[,MLE]*sum(ODD@gmax$gmax)/mean(sum(Disp[,MLE])) %>% round()
  ODD@data$Disp<-Dam[,MLE,1]  #should this be named data$DispPred or something instead?
  ODD@data$Mort<-Dam[,MLE,2]
  ODD@data$nBD<-Dam[,MLE,3]
  # I know this is kind of repeating code, but I want the other function as fast as possible
  ODD@predictDisp<-funcy(MLE,LLout=F) 
  
  return(ODD)
  
})
