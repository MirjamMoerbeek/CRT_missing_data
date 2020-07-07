source("functions.R")

server <- function(input, output) {

###########################################################################################################################################################################
### attrition graphs
###########################################################################################################################################################################

  output$survivalplot <- renderPlot({
    
    omega0=input$omega0
    gamma0=input$gamma0
    maxdur=input$maxweeks*7
    time=seq(0,1,length=maxdur)
    remain0 = (1 - omega0)^time^gamma0
    
    omega1=input$omega1
    gamma1=input$gamma1
    remain1 = (1 - omega1)^time^gamma1
    par(mfrow=c(1,2))
    plot(seq(1,maxdur),remain0,type="b",pch=16,ylim=c(0,1),xlim=c(0,maxdur),xlab="time in days",ylab="survival probability function",main="Control condition")
    plot(seq(1,maxdur),remain1,type="b",pch=16,ylim=c(0,1),xlim=c(0,maxdur),xlab="time in days",ylab="survival probability function",main="Intervention condition")
  })
  
  output$hazardplot <- renderPlot({
    omega0=input$omega0
    gamma0=input$gamma0
    maxdur=input$maxweeks*7
    time=seq(0,1,length=maxdur)
    remain0 = (1 - omega0)^time^gamma0
    hazard0=rep(0,maxdur-1)
    for(kk in 1:maxdur-1)    
      hazard0[kk]=(remain0[kk]-remain0[kk+1])/remain0[kk]
    if(gamma0==1) hazard0=rep(hazard0[1],maxdur-1)

    omega1=input$omega1
    gamma1=input$gamma1
    remain1 = (1 - omega1)^time^gamma1
    hazard1=rep(0,maxdur-1)
    for(kk in 1:maxdur-1)    
      hazard1[kk]=(remain1[kk]-remain1[kk+1])/remain1[kk]
    if(gamma1==1) hazard1=rep(hazard1[1],maxdur-1)

    max.hazard=max(c(hazard0,hazard1))
    par(mfrow=c(1,2)) 
    plot(seq(1,maxdur-1),hazard0,type="b",pch=16,ylim=c(0,max.hazard),xlim=c(0,maxdur),xlab="time in days",ylab="hazard probability function",main="Control condition")
    plot(seq(1,maxdur-1),hazard1,type="b",pch=16,ylim=c(0,max.hazard),xlim=c(0,maxdur),xlab="time in days",ylab="hazard probability function",main="Intervention condition")
  })

###########################################################################################################################################################################
### Results graphs
###########################################################################################################################################################################
  
  output$Resultsplot <- renderPlot({
    
    validate(need(input$maxweeks>=input$R1,"Number of weeks in Design 1 cannot exceed the maximum duration of the trial")
    )
    validate(need(input$maxweeks>=input$R2,"Number of weeks in Design 2 cannot exceed the maximum duration of the trial")
    )
    validate(need(input$maxweeks>=input$R3,"Number of weeks in Design 3 cannot exceed the maximum duration of the trial")
    )
    validate(need(input$maxweeks>=input$R4,"Number of weeks in Design 4 cannot exceed the maximum duration of the trial")
    )
    validate(need(input$maxweeks>=input$R5,"Number of weeks in Design 5 cannot exceed the maximum duration of the trial")
    )
    
    m.min=input$m[1]
    m.max=input$m[2]
    
    decay=input$decay
    ICC=input$ICC
    
    delta=input$delta
    alpha=input$alpha
    test=input$test
    if(test=="2") alpha=alpha/2

    maxdays=input$maxweeks*7
    omega0=input$omega0
    gamma0=input$gamma0
    omega1=input$omega1
    gamma1=input$gamma1
    
    k1=input$k1
    R1=input$R1
    weekdays1=as.numeric(input$checkDays1)
    
    k2=input$k2
    R2=input$R2
    weekdays2=as.numeric(input$checkDays2)
    
    k3=input$k3
    R3=input$R3
    weekdays3=as.numeric(input$checkDays3)
    
    k4=input$k4
    R4=input$R4
    weekdays4=as.numeric(input$checkDays4)
    
    k5=input$k5
    R5=input$R5
    weekdays5=as.numeric(input$checkDays5)

    m.vec=seq(m.min,m.max)
    var1.vec=seq(m.min,m.max)
    var2.vec=seq(m.min,m.max)
    var3.vec=seq(m.min,m.max)
    var4.vec=seq(m.min,m.max)
    var5.vec=seq(m.min,m.max)

    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', value = 0, {
    
    for(ii in 1:length(m.vec))
    {
      var1.vec[ii]=f.varCRT(r=1-decay,rho=ICC,m=m.vec[ii],k=k1,R=R1,maxdays=maxdays,weekdays=weekdays1,omega0=omega0,gamma0=gamma0,omega1=omega1,gamma1=gamma1)
      var2.vec[ii]=f.varCRT(r=1-decay,rho=ICC,m=m.vec[ii],k=k2,R=R2,maxdays=maxdays,weekdays=weekdays2,omega0=omega0,gamma0=gamma0,omega1=omega1,gamma1=gamma1)
      var3.vec[ii]=f.varCRT(r=1-decay,rho=ICC,m=m.vec[ii],k=k3,R=R3,maxdays=maxdays,weekdays=weekdays3,omega0=omega0,gamma0=gamma0,omega1=omega1,gamma1=gamma1)
      var4.vec[ii]=f.varCRT(r=1-decay,rho=ICC,m=m.vec[ii],k=k4,R=R4,maxdays=maxdays,weekdays=weekdays4,omega0=omega0,gamma0=gamma0,omega1=omega1,gamma1=gamma1)
      var5.vec[ii]=f.varCRT(r=1-decay,rho=ICC,m=m.vec[ii],k=k5,R=R5,maxdays=maxdays,weekdays=weekdays5,omega0=omega0,gamma0=gamma0,omega1=omega1,gamma1=gamma1)
      incProgress(1/length(m.vec))
      Sys.sleep(0.05)
    }
                 })
    
    max.var=max(c(var1.vec,var2.vec,var3.vec,var4.vec,var5.vec))
    max.RE=max(var1.vec/var2.vec,var1.vec/var3.vec,var1.vec/var4.vec,var1.vec/var5.vec)
    if(max.RE<1) max.RE=1
    
    power1.vec=pnorm(delta/sqrt(var1.vec)-qnorm(1-alpha/2))
    power2.vec=pnorm(delta/sqrt(var2.vec)-qnorm(1-alpha/2))
    power3.vec=pnorm(delta/sqrt(var3.vec)-qnorm(1-alpha/2))
    power4.vec=pnorm(delta/sqrt(var4.vec)-qnorm(1-alpha/2))
    power5.vec=pnorm(delta/sqrt(var5.vec)-qnorm(1-alpha/2))
    par(mfrow=c(1,3)) 
   
   ### Variance graphs
   plot(m.vec,var1.vec,type="l",col="black",lty=2,xlim=c(0,m.max),ylim=c(0,max.var),xlab="number of subjects per cluster-period",ylab="variance",main="Variance of the treatment effect estimator",cex.axis=2,cex.lab=1.75,cex.main=1.75,lwd=2)
   lines(m.vec,var2.vec,col="red",lwd=2,lty=3)
   lines(m.vec,var3.vec,col="blue",lwd=2,lty=1)
   lines(m.vec,var4.vec,col="darkorange",lwd=2,lty=4)
   lines(m.vec,var5.vec,col="green",lwd=2,lty=5)
   legend("topright", legend=c("Design 1","Design 2","Design 3","Design 4","Design 5"),lty=c(2,3,1,4,5),col=c("black","red","blue","darkorange","green"),box.lty=0,cex=1.75,bg="transparent",lwd=2)
   
   ### Power graphs
   plot(m.vec,power1.vec,type="l",col="black",lty=2,xlim=c(0,m.max),ylim=c(0,1),xlab="number of subjects per cluster-period",ylab="power",main="Power of the test on treatment effect", cex.axis=1.75,cex.lab=1.75,cex.main=1.75,lwd=2)
   lines(m.vec,power2.vec,col="red",lwd=2,lty=3)
   lines(m.vec,power3.vec,col="blue",lwd=2,lty=1)
   lines(m.vec,power4.vec,col="darkorange",lwd=2,lty=4)
   lines(m.vec,power5.vec,col="green",lwd=2,lty=5)
   legend("bottomright", legend=c("Design 1","Design 2","Design 3","Design 4","Design 5"),lty=c(2,3,1,4,5),col=c("black","red","blue","darkorange","green"),box.lty=0,cex=1.75,bg="transparent",lwd=2)
    
    ### Eficiency relative to design 1
    plot(m.vec,var1.vec/var1.vec,type="l",col="black",lty=2,xlim=c(0,m.max),ylim=c(0,max.RE),xlab="number of subjects per cluster-period",ylab="efficiency relative to Design 1",main="Relative Efficiency", cex.axis=1.75,cex.lab=1.75,cex.main=1.75,lwd=2)
    lines(m.vec,var1.vec/var2.vec,col="red",lwd=2,lty=3)
    lines(m.vec,var1.vec/var3.vec,col="blue",lwd=2,lty=1)
    lines(m.vec,var1.vec/var4.vec,col="darkorange",lwd=2,lty=4)
    lines(m.vec,var1.vec/var5.vec,col="green",lwd=2,lty=5)
    legend("bottomright", legend=c("Design 1","Design 2","Design 3","Design 4","Design 5"),lty=c(2,3,1,4,5),col=c("black","red","blue","darkorange","green"),box.lty=0,cex=1.75,bg="transparent",lwd=2)
  })
  
  
  ###########################################################################################################################################################################
  ### Results table
  ###########################################################################################################################################################################
  
  output$ResultsTable <- renderDataTable({
    
    validate(need(input$maxweeks>=input$R1,"Number of weeks in Design 1 cannot exceed the maximum duration of the trial")
    )
    validate(need(input$maxweeks>=input$R2,"Number of weeks in Design 2 cannot exceed the maximum duration of the trial")
    )
    validate(need(input$maxweeks>=input$R3,"Number of weeks in Design 3 cannot exceed the maximum duration of the trial")
    )
    validate(need(input$maxweeks>=input$R4,"Number of weeks in Design 4 cannot exceed the maximum duration of the trial")
    )
    validate(need(input$maxweeks>=input$R5,"Number of weeks in Design 5 cannot exceed the maximum duration of the trial")
    )
    
    n=input$nrclusters
    m=input$nrsubjects
    p=input$nrperiods
    
    corr=input$basecorr
    r=1-input$decay
    
    c=input$costcluster
    s=input$costsubject
    x=input$costcrossover
    budget=input$budget

    m.min=input$m[1]
    m.max=input$m[2]
    
    decay=input$decay
    ICC=input$ICC

    delta=input$delta
    alpha=input$alpha
    test=input$test
    if(test=="2") alpha=alpha/2

    maxdays=input$maxweeks*7
    omega0=input$omega0
    gamma0=input$gamma0
    omega1=input$omega1
    gamma1=input$gamma1
    
    k1=input$k1
    R1=input$R1
    weekdays1=as.numeric(input$checkDays1)
    
    k2=input$k2
    R2=input$R2
    weekdays2=as.numeric(input$checkDays2)
    
    k3=input$k3
    R3=input$R3
    weekdays3=as.numeric(input$checkDays3)
    
    k4=input$k4
    R4=input$R4
    weekdays4=as.numeric(input$checkDays4)
    
    k5=input$k5
    R5=input$R5
    weekdays5=as.numeric(input$checkDays5)
    
    m.vec=seq(m.min,m.max)
    var1.vec=seq(m.min,m.max)
    var2.vec=seq(m.min,m.max)
    var3.vec=seq(m.min,m.max)
    var4.vec=seq(m.min,m.max)
    var5.vec=seq(m.min,m.max)
    
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', value = 0, {
                   
                   for(ii in 1:length(m.vec))
                   {
                     var1.vec[ii]=f.varCRT(r=1-decay,rho=ICC,m=m.vec[ii],k=k1,R=R1,maxdays=maxdays,weekdays=weekdays1,omega0=omega0,gamma0=gamma0,omega1=omega1,gamma1=gamma1)
                     var2.vec[ii]=f.varCRT(r=1-decay,rho=ICC,m=m.vec[ii],k=k2,R=R2,maxdays=maxdays,weekdays=weekdays2,omega0=omega0,gamma0=gamma0,omega1=omega1,gamma1=gamma1)
                     var3.vec[ii]=f.varCRT(r=1-decay,rho=ICC,m=m.vec[ii],k=k3,R=R3,maxdays=maxdays,weekdays=weekdays3,omega0=omega0,gamma0=gamma0,omega1=omega1,gamma1=gamma1)
                     var4.vec[ii]=f.varCRT(r=1-decay,rho=ICC,m=m.vec[ii],k=k4,R=R4,maxdays=maxdays,weekdays=weekdays4,omega0=omega0,gamma0=gamma0,omega1=omega1,gamma1=gamma1)
                     var5.vec[ii]=f.varCRT(r=1-decay,rho=ICC,m=m.vec[ii],k=k5,R=R5,maxdays=maxdays,weekdays=weekdays5,omega0=omega0,gamma0=gamma0,omega1=omega1,gamma1=gamma1)
                     incProgress(1/length(m.vec))
                     Sys.sleep(0.05)
                   }
                 })
    
    power1.vec=pnorm(delta/sqrt(var1.vec)-qnorm(1-alpha/2))
    power2.vec=pnorm(delta/sqrt(var2.vec)-qnorm(1-alpha/2))
    power3.vec=pnorm(delta/sqrt(var3.vec)-qnorm(1-alpha/2))
    power4.vec=pnorm(delta/sqrt(var4.vec)-qnorm(1-alpha/2))
    power5.vec=pnorm(delta/sqrt(var5.vec)-qnorm(1-alpha/2))
    
    RE1.vec=var1.vec/var1.vec
    RE2.vec=var1.vec/var2.vec
    RE3.vec=var1.vec/var3.vec
    RE4.vec=var1.vec/var4.vec
    RE5.vec=var1.vec/var5.vec

    data=round(cbind(seq(m.min,m.max),var1.vec,var2.vec,var3.vec,var4.vec,var5.vec,power1.vec,power2.vec,power3.vec,power4.vec,power5.vec,RE1.vec,RE2.vec,RE3.vec,RE4.vec,RE5.vec),5)
    data=as.data.frame(data)
    colnames(data)<-c("m","Variance 1","Variance 2","Variance 3","Variance 4","Variance 5","Power 1", "Power 2", "Power 3","Power 4","Power 5","RE1 ","RE 2","RE 3","RE 4","RE 5")
    data
  })
  
}