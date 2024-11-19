qbeta(0.025,5,3)
qbeta(0.975,5,3)
1-pbeta(0.5,5,3)

a<-b<-0
awins<-1-(gamma(a+b+8)/(gamma(a+5)*gamma(b+3)))/(gamma(a+b+11)/(gamma(a+5)*gamma(b+6)))

awins/(1-awins)

1/18
1-1/11

485/512

5/8

pdf("/Users/sib2/Box Sync/Faculty/Education/BIOSTAT725/website/slides/images/01/beta64.pdf", width = 5, height = 5)
x<-seq(0,1,0.001)
density1<-dbeta(x,6,4)
plot(x,density1,type="l",ylim=c(0,3),ylab="Density",xlab=expression(pi),col="red", lwd = 1.5)
legend("topleft",legend=c("Beta(6,4)"),lty=c(1),col=c("red"), lwd = 1.5, bty = "n")
abline(v=qbeta(0.025,6,4), lwd = 1.5)
abline(v=qbeta(0.975,6,4), lwd = 1.5)
dev.off()


pdf("/Users/Sam/Documents/Sam/School/Graduate/DCRI/Talk09162015/beta.pdf")
x<-seq(0,1,0.001)
density1<-dbeta(x,1,1)
density0<-x^(-1)*(1-x)^(-1)
plot(x,density0/10,type="l",ylim=c(0,10),ylab="Density",xlab=expression(pi),col="red")
lines(x,density1,type="l",col="green")
legend(0.5,10,legend=c("Beta(1,1)","Beta(0,0)"),lty=c(1,1),col=c("green","red"))
dev.off()

awins<-1-(3/8)*(4/9)*(5/10)

awins/(1-awins)