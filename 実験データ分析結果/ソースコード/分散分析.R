##### 分散分析ＡｓＢＣ(３要因混合)：主分析

# パッケージ car   が必要
# パッケージ psych が必要

## js-STARの入力: AsBC
levA = 2          # level of A
levB = 2          # level of B
levC = 8          # level of C

n    = c( # データ数
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8
)

data = c(
  678.4314444, 659.1594, 686.340875, 776.249875, 940.5086364, 830.376875, 838.8617, 674.0702, 664.1602, 718.6206667, 756.4770909, 834.9097273, 888.4891, 906.436, 789.0171111, 736.989125,
  1016.495, 1024.955111, 934.8521, 954.813625, 1256.737333, 969.5596667, 1037.295455, 1001.6117, 955.1493, 1025.0048, 951.6938182, 1342.935273, 1117.881, 1305.2297, 1057.5783, 1117.5105,
  1349.8328, 1495.578182, 2032.616333, 1975.171222, 2196.188333, 1966.564, 1553.567455, 1581.7711, 1076.8674, 1670.502889, 1770.0383, 2013.434333, 2165.102222, 1867.695545, 1853.517667, 1662.569727,
  969.1294444, 1071.0668, 1598.621667, 1448.57, 1676.807167, 2039.1725, 1262.742667, 1060.993375, 1434.826625, 1283.326667, 1583.289111, 1390.2165, 1853.169714, 1568.7612, 1586.431714, 1177.581571,
  1025.1724, 915.6491998, 1031.703, 1253.9664, 1703.003428, 1484.855857, 1053.423667, 969.0243998, 868.8217997, 943.1025998, 1215.617, 1435.020555, 2158.304667, 1418.4654, 1065.8986, 939.5019998,
  745.2395, 879.4888182, 850.8387, 1163.121222, 1267.761909, 1195.362444, 897.5007273, 738.3202, 761.2539, 937.4703, 949.551, 1574.965444, 1308.025625, 1544.907875, 1023.647, 847.3760909,
  610.7886, 735.9464545, 983.964, 951.916, 1145.473143, 830.614, 836.2018182, 643.7025, 810.4569, 1138.7325, 990.2003636, 1198.272444, 1387.951, 991.5297273, 859.8574, 875.2741,
  541.5, 530.1818182, 654.7, 663.5, 938.5714286, 945.25, 675.2, 605.7, 668, 712.5, 853.5454545, 898.8181818, 894.5, 851.2727273, 622.7, 726.9090909,

  1062.1278, 1279.8897, 998.2982, 1315.5856, 1319.199, 1290.0609, 1029.289444, 1166.6383, 1313.5543, 1606.686333, 1434.804222, 1882.699444, 1966.026, 1682.976571, 1579.5035, 1239.969222,
  913.6425, 932.3578889, 976.5072222, 1023.906375, 1064.126286, 1097.488714, 940.364, 897.1232, 1057.022667, 1012.552778, 1146.726, 1084.573667, 1160.450556, 1179.261875, 1078.8407, 1189.533375,
  480.2758, 459.5651, 519.7992222, 651.9268, 878.9955556, 732.734875, 571.5995556, 493.2409, 476.8002, 514.6099, 557.7712, 695.7194, 848.0311111, 742.1257, 595.9416667, 551.6991,
  781.84725, 648.39925, 820.1438889, 867.514875, 950.05, 912.703, 903.777625, 744.6812857, 1104.77975, 949.2632222, 946.7183333, 1031.042143, 937.5535, 912.211125, 1104.735429, 808.574375,
  589.9158, 563.0528889, 618.8024, 711.7395, 760.638, 910.2153333, 724.6714444, 589.323, 549.9149, 654.5257, 700.046, 792.3297, 901.5985, 808.788, 634.2434, 672.4436,
  781.2098889, 847.6823, 683.592, 862.6418889, 1071.736778, 989.3222, 736.5058, 840.2935555, 1003.1178, 824.959, 900.8591, 1071.199625, 1892.9565, 990.160125, 870.8914444, 908.8337778,
  984.8149, 717.6304, 753.6392222, 1152.520667, 1281.181889, 930.833, 769.651, 718.9001, 706.9624444, 838.7693, 2050.824, 1025.0947, 1094.572, 956.7593333, 902.8592, 808.4275,
  695.867, 708.5347, 787.9246, 929.6885, 1147.527143, 1084.644125, 775.4813, 722.8359, 791.7885, 830.6436, 851.3037, 1097.2495, 1191.833889, 909.5507, 923.3105, 844.1068889
)

taju="holm"  # BH, holm, b

## スタック
A=c(); B=c(); s=c(); C=c()
jun=0; gun=0
for(i in 1:levA){ gun=gun+1
for(j in 1:levB){ for(k in 1:levC){ jun=jun+1
for(L in 1:n[jun]) A=c(A,i)
} } }
A=factor(A)
N= sum(n)/levB/levC
s= gl(N, levB*levC, sum(n) )
B= gl(levB, levC, sum(n) )
C= gl(levC, 1, sum(n) )
AsBC= data.frame(A, s, B, C, data )

## 基本統計
Title="ＡxＢxＣの平均とSD(不偏分散の平方根)"
JUNB= function(L) if(L<4) sum(1:L) else floor(L*(1+sign(L-4)*0.5))
nmat= matrix(n, nc=levB*levC, by=1 )
library( car )
library( psych )
tx0=c(); txD=c(); mds=c()
for( i in 1:levA ){ tx0=c(tx0,rep(NA,5)); txD=rbind(txD,rep(NA,5)); mds=c(mds,"")
for( j in 1:levB ){ for( k in 1:levC ){
  kx= describe(subset(AsBC,A==i & B==j & C==k)$data )
  tx0= c(tx0, kx$n,kx$mean,kx$sd,kx$min,kx$max )
  mds= c(mds, paste("A.",i,"_B.",j,"_C.",k, sep="") )
  txD= rbind(txD, as.numeric(kx[1,c(5,10:13)]) )
}  }  }
tx0= matrix(round(tx0,4), nc=5,by=1 )
tx0= tx0[-1,]
rownames(tx0)= mds[-1]
colnames(tx0)= c("   n","    Mean","     SD","   Min","   Max" )
txD= round(txD[-1,], 2)
rownames(txD)= rownames(tx0)
colnames(txD)= c("Median","Range","歪度","尖度","SE" )
hk= tx0[,2][ is.na(tx0[,2])==0 ]
sd= tx0[,3][ is.na(tx0[,3])==0 ]

# 作図
mitsdo= c(0:(levC-1))*10*log(1:levC) # バー斜線密度
mitsdo= rep( mitsdo, levA*levB )     # Ｃ水準を林立
kankak= c()
for(i in 1:levA){ for(j in 1:levB){
  if(j==1) spc=2 else spc=1
  kankak= c(kankak, spc, rep(0,levC-1))
}  }
bai = 2
shoh= min(min(hk),0) # 最小平均か，ゼロ
dais= max(sd)        # 最大ＳＤ
yjik= c(shoh+dais*bai*sign(shoh), max(hk)+dais*bai ) # ｙ軸範囲

for( L in 1:2 )
{
  if( L==1 ) Sd=sd else { Sd=txD[,5][is.na(txD[,5])==0]; Title="ＡxＢxＣの平均とSE(標準誤差)" }
  quartz()
  antx= barplot(hk, yli=yjik, # ｘ軸の範囲指定しない
                las=1, tck=0.02,  # ｙ軸の目盛デザイン
                col=1, den=mitsdo,
                spa=kankak,       # バーの間隔
                names.arg=c(rep(paste("C",1:levC,sep=""),levA*levB)),
                xlab="", ylab="平　均", # ラベル名
                sub=paste("Ｂ1(左端",levC,"本で１水準) ～Ｂ",levB,"　in　Ａ1(左端から",levC*levB,"本で１水準) ～Ａ",levA, sep="" ),
                main=Title,
                cex.main=1.5, cex.sub=1.2, cex.lab=1.2 )
  arrows(antx,hk,antx,hk+Sd*sign(hk),lwd=L,ang=90,len=0.12 )
  if(L==2) arrows(antx,hk,antx,hk-Sd*sign(hk),lwd=L,ang=90,len=0.12 )
}

# Type3
library( car )
TsF= function(ss,d1,d2,MS,adj=0)
{
  ms=ss/d1; fc=ms/MS
  pc= pf(fc,d1,d2,low=0)
  if(adj==1) pc=p.adjust(pc,me=taju)
  et= fc*d1/(fc*d1+d2)
  matrix(round(c(ss,d1,ms,fc,pc,et),4),nc=6 )
}
NAM= c("主効果Ａ","主効果Ｂ","主効果Ｃ","ＡxＢ","ＡxＣ","ＢxＣ","二次の交互作用ＡxＢxＣ" )

ty1= anova(lm(data~A*B*C+s+s:B+s:C,AsBC) )
kx3= lm(data~A*B*C,AsBC,contrasts=list(A=contr.sum,B=contr.sum,C=contr.sum) )
ty3= Anova(kx3,type=3 )

DFa=levA-1;DFb=levB-1;DFc=levC-1;DFs=N-levA
SS= c(ty3$S[2],ty1$S[4],ty3$S[c(3,5)],ty1$S[8],ty3$S[c(4,6)],ty1$S[9],ty3$S[c(7,8)],ty1$S[11] )
DF= c(DFa,DFs, DFb,DFa*DFb,DFs*DFb, DFc,DFa*DFc,DFs*DFc, DFb*DFc,DFa*DFb*DFc,DFs*DFb*DFc )
d2= c(DF[2],NA,DF[c(5,5)],NA,DF[c(8,8)],NA,DF[c(11,11)],NA )
MS= SS/DF
ms= c(MS[2],NA,MS[c(5,5)],NA,MS[c(8,8)],NA,MS[c(11,11)],NA )

tx1= c()
tx1= TsF(SS,DF,d2,ms)
Fch=tx1[,4]; Pch=tx1[,5]; eta=tx1[,6]
rownames(tx1)= c(" 主効果Ａ","　　ｓ"," 主効果Ｂ","　Ａ×Ｂ","　ｓｘＢ"," 主効果Ｃ","　Ａ×Ｃ","　ｓｘＣ","　Ｂ×Ｃ","Ａ×Ｂ×Ｃ","ｓｘＢｘＣ" )
colnames(tx1)= c("TypeⅢ_SS ","df ", "MS ", "Ｆ ", "ｐ ","  ηp2 " )

kx= Anova(kx3,type=2 ) # Type2
ss= c(kx$S[1],SS[2],kx$S[c(2,4)],SS[5],kx$S[c(3,5)],SS[8],kx$S[c(6,7)],SS[11] )
ty2= c()
ty2= TsF(ss,DF,d2,ms)
rownames(ty2)= rownames(tx1)
colnames(ty2)= c("TypeⅡ_SS ","df ", "MS ", "Ｆ ", "ｐ ","  ηp2 " )

# 分散の均一性
tx8=c();mds=c(); Bst=c();Bdf=c();Bpv=c();Bsj=c()
for(j in 1:levB){ for(k in 1:levC){
  kx8= bartlett.test(data~A,subset(AsBC,B==j & C==k) )
  Bst=c(Bst,kx8$stat); Bdf=c(Bdf,kx8$para)
  pch=kx8$p.va; Bpv=c(Bpv,pch)
  if( pch<0.05 ) Bsj=c(Bsj, j, k) # 有意の水準番号
  mds= c(mds, paste("at_B",j,"_C",k, sep="") )
}  }
tx8= matrix(round(c(Bst,Bdf,Bpv),4), nc=3 )
rownames(tx8)= mds
colnames(tx8)= c("    χ2"," df","     ｐ" )

# 球面性
HSOK= function( Lev )
{
  rmat= cor(dmat)
  keis= c()
  for(x in 1:(Lev-1)){ for(y in (x+1):Lev) keis=c(keis,rmat[x,y] ) }
  tanh(sum(nrow(dmat)*atanh(keis))/( nrow(dmat)*length(keis) ) ) # Fisher.Z
}
MIXR= function(rmat,lev) if(sum(rmat<0)==0 | (sum(rmat<0)>0 & sum(rmat>0)==lev) ) 0 else 1

dmat=c(); mixR=c()
for(j in 1:levB) dmat=cbind(dmat,subset(AsBC$data,B==j) )
corB= HSOK( levB )
Bmat=cor(dmat);mds=paste("B",1:levB,sep="");dimnames(Bmat)=list(mds,mds)
mixR= c(mixR,MIXR(Bmat,levB) )

dmat= c()
for(k in 1:levC) dmat=cbind(dmat,subset(AsBC$data,C==k) )
corC= HSOK( levC )
Cmat=cor(dmat);mds=paste("C",1:levC,sep="");dimnames(Cmat)=list(mds,mds)
mixR= c(mixR,MIXR(Cmat,levC) )

dmat= c()
for(j in 1:levB){ for(k in 1:levC) dmat=cbind(dmat,subset(AsBC$data,B==j & C==k)) }
corT= HSOK( levB*levC )
Tmat=cor(dmat);mds= paste("B",rep(1:levB,each=levC),"_C",rep(1:levC, levB),sep="")
dimnames(Tmat)= list(mds, mds )
mixR= c(mixR,MIXR(Tmat,levB*levC) )

AA=c()
for(i in 1:levA) AA=c(AA,rep(i,nmat[i,1]) )
AA= factor(AA)
lmod= lm(dmat~AA)  # 線形モデル
idat= data.frame(sC=gl(levC,1,levB*levC),sB=gl(levB,levC) ) # 逆階層

kxB= mauchly.test(lmod, M=~sC+sB,X=~sC, idata=idat )
kxC= mauchly.test(lmod, M=~sC+sB,X=~sB, idata=idat )
kxX= mauchly.test(lmod, X=~sC+sB, idata=idat )
MsW= c(kxB$stat, kxC$stat, kxX$stat )
MsP= c(kxB$p.va, kxC$p.va, kxX$p.va )
epGG=c(); epHF=c(); goku="Spherical"
eps= c(attributes(anova(lmod,M=~sB+sC,X=~sC,test=goku,idata=idat))$head[8:9],
       attributes(anova(lmod,M=~sB+sC,X=~sB,test=goku,idata=idat))$head[8:9],
       attributes(anova(lmod,      X=~sB+sC,test=goku,idata=idat))$head[5:6] )
eps= as.numeric(substr(eps,29,34) )
epGG= eps[c(1,3,5)]
epHF= eps[c(2,4,6)]

tx9= c()
tx9= matrix(round(c(MsW,MsP,epGG,epHF), 4), nr=3 )
rownames(tx9)= c("要因Ｂ","要因Ｃ","Ｂ×Ｃ" )
colnames(tx9)= c("Mauchly's W","   ｐ値","  G-G_ε"," H-F_ε" )

crP= c()
epGG=rep(epGG,each=2); epHF[ epHF>1 ]=1; epHF=rep(epHF,each=2)
gyo= c(3,4, 6,7, 9,10 )
GGp= pf(Fch[gyo],DF[gyo]*epGG,d2[gyo]*epGG,low=0 )
HFp= pf(Fch[gyo],DF[gyo]*epHF,d2[gyo]*epHF,low=0 )
crP= matrix(round(c(DF[gyo],Fch[gyo],Pch[gyo],GGp,HFp),4),nc=5 )
rownames(crP)= c("主効果Ｂ"," ＡxＢ","主効果Ｃ"," ＡxＣ"," ＢxＣ","ＡxＢxＣ" )
colnames(crP)= c("df","     Ｆ","     ｐ"," G-G_ｐ"," H-F_ｐ" )

# パワ
gyo= c(1,gyo)
ESf= sqrt(Fch[gyo]*DF[gyo]/d2[gyo] )

NCPz= c( # 相関＝zero
  ESf[1]^2*N*(DF[9]+1),
  ESf[2:3]^2*N*levB,
  ESf[4:5]^2*N*levC,
  ESf[6:7]^2*N*((levB-1)*(levC-1)+1)
)
tmp= c(corT,corB,corB,corC,corC,corT,corT )
corT=abs(corT);corB=abs(corB);corC=abs(corC)
NCPr= c(
  ESf[1]^2*N*(DF[9]+1)*( 1/(DF[9]*corT+(1-corT)) ),
  ESf[2:3]^2*N*(levB)*( 1/(1-corB) ),
  ESf[4:5]^2*N*(levC)*( 1/(1-corC) ),
  ESf[6:7]^2*N*( (levB-1)*(levC-1)+1 )*( 1/(1-corT) )
)
crF= qf(0.95,DF[gyo],d2[gyo] )
PWz= 1-pf(crF,DF[gyo],d2[gyo],NCPz )
PWr= 1-pf(crF,DF[gyo],d2[gyo],NCPr )

tx2= c()
tx2= matrix(round(c(ESf,PWz,PWr,tmp),4), nr=7 )
rownames(tx2)= c(" 主効果Ａ"," 主効果Ｂ","　Ａ×Ｂ"," 主効果Ｃ","　Ａ×Ｃ","　Ｂ×Ｃ","Ａ×Ｂ×Ｃ" )
colnames(tx2)= c("効果量ｆ"," 検出力0"," 検出力r"," 水準間相関")

# 主効果
hkA=c();hkB=c();hkC=c();mcA=c();mcB=c();mcC=c()
Xari= matrix(c(4,7,10, 4,9,10, 7,9,10), nc=3, by=1 )
levs= c(levA,levB,levC)
nama= c("A","B","C" )

kx= anova(lm(AsBC[,5]~AsBC$A ))
plSD= matrix(c(sqrt(kx$M[2]),kx$D[2]), nc=2 )
rownames(plSD)= NAM[1]
colnames(plSD)= c("pooled_SD","   DF" )

for( L in 1:3 )
{
  hyo1=c(); hyo2=c()
  x=c(); y=c(); u=c()
  for( i in 1:levs[L] )
  {
    dt= subset(AsBC$data,AsBC[,round(L*1.4)]==i )
    u=c(u, length(dt) )
    x=c(x, mean(dt) )
    y=c(y, sd(dt) )
  }
  hyo1= rbind(matrix(round(c(u,x,y),4),nr=3,by=1),rep(NA,levs[L]) )
  rownames(hyo1)= c("ｎ","Mean","S.D.","-----" )
  colnames(hyo1)= paste(nama[L],1:levs[L], sep="" )
  if( levs[L]>2 )
  {
    hyo2= round(pairwise.t.test(AsBC[,5],AsBC[,sign(L-1)+L],p.ad=taju,pair=sign(L-1))$p.v, 4 )
    rownames(hyo2)= paste(nama[L],2:levs[L], sep="" )
    colnames(hyo2)= paste(nama[L],1:(levs[L]-1), sep="" )
  }
  if(L==1){ hkA=hyo1; mcA=hyo2 }
  if(L==2){ hkB=hyo1; mcB=hyo2 }
  if(L==3){ hkC=hyo1; mcC=hyo2 }
}


# ■結果の書き方
PV= function(x) floor(x*1000)/1000
txt=c()
if( sum(sd==0)>0 ) txt="※分散=0 の群があります。全データが同じ値です。以下の分析は信頼性がありません。\n\n"
txt= paste(txt, "　各群の各水準の○○得点について基本統計量をTable(tx0)に示す。\n　要因Ａを参加者間，要因Ｂ・Ｃを参加者内に配置した３要因分散分析（TypeⅢ_SS使用）を行った結果 (Table(tx1)参照)，", sep="" )

Pow= c()
for(gyo in 1:7)
{
  if( mixR[(round(1/gyo)*6+gyo)%/%2]<1 ) Pow=c(Pow,PWr[gyo]) else Pow=c(Pow,PWz[gyo])
}
Pow= Pow[c(1,2,4,3,5:7)]

yui=c()
for( L in 1:7 )
{
  if( L==4 ) txt=paste(txt, "また一次の交互作用については，", sep="" )
  JB= JUNB(L)
  pch= Pch[JB]

  if( pch<0.05 ){ yui=c(yui,L); if(L==3 | L>5 ) tmp="であった" else tmp="であり" } else
    if( pch<0.10 ){ yui=c(yui,L); if(L==3 | L>5 ) tmp="傾向であった" else tmp="傾向であり" } else
    { if(L==3 | L>5 ) tmp="でなかった" else tmp="でなく" }
  if(L==3 | L>5 ) goku="。" else goku="，"
  fuki= paste(" (F(",DF[JB],",",d2[JB],")=",round(Fch[JB],3),", p=",PV(pch),", ηp2=",round(eta[JB],3),", 1-β=",round(Pow[L],3),")",goku, sep="" )
  if(L==7) goku="は" else goku="が"
  if(L==7 & pch<0.10 & length(yui)>1 ){ mata="ただし，"; goku="が" } else mata=""
  txt = paste(txt, mata,NAM[L],goku,"有意",tmp,fuki, sep="" )
}

# 検出力
kaz=length(yui); yosh=0; mae=9
if( kaz>0 ){ for( i in 1:kaz ){
  pwr= round(Pow[yui[i]], 3)
  if( pwr>0.80 ){ tmp="十分である。"; yosh=1 } else
    if( pwr>0.76 ){ tmp="ほぼ十分といえる。"; yosh=1 } else
      if( pwr>0.70 ){ tmp="やや低いが0.70以上あり不十分ではない。"; yosh=0 } else { tmp="不十分であり，信頼性が低い。"; yosh=0 }
  if( mae==yosh ) goku="も" else goku="は"
  mae= yosh
  if( i==1 ){ mata="\n　"; fuki=" (1-β) " } else { mata=""; fuki="" }
  txt= paste(txt, mata,NAM[yui[i]],"の検出力",fuki,goku,tmp, sep="" )
}  }
tmp="正負が混在している場合は平均相関を0と仮定し，それ以外はFisherの重み付きZ変換値による平均相関を"
if(levB>2 | levC>2) goku="に" else { goku="を"; tmp="" }
txt= paste(txt, "なお検出力の値は水準間の相関係数",goku,tmp,"用いて算出した。", sep="" )

# 分散均一性
if( min(Pch, na.rm=1)<0.10 )
{
  txt= paste(txt, "\n　参加者間要因の分散の均一性についてBartlett検定を行った結果 (Table(tx8)参照)，要因Ｂ・Ｃの", sep="" )
  kaz=0; gok1="s<"; gok2="s>"
  pch=min(Bpv); kai=round(max(Bst),3); goku=""
  if( pch>0.05 ) tmp="いずれの水準においても有意でないことを確認した" else
  {
    kaz= length(Bsj)/2
    tmp=c(); gok0=""; sjn=c()
    for(z in 1:kaz){ if(z>1) gok0="，"; j=(z-1)*2+1; k=(z-1)*2+2
    tmp=paste(tmp, gok0,"B",Bsj[j],"_C",Bsj[k], sep=""); sjn=c(sjn, (Bsj[j]-1)*levC+Bsj[k] ) }
    tmp= paste("水準",tmp,"において有意であった", sep="" )
    if( kaz==1 ){ gok1="="; gok2="=" }
    if( kaz>1 ){ gok1="s>"; gok2="s<"; pch=PV(max(Bpv[sjn])); kai=round(min(Bst[sjn]),3) }
    goku="以下，参考までに分析を進める。"
  }
  fuki= paste(" (χ2(",Bdf[1],")",gok1,kai,", p",gok2,PV(pch),")。", sep="" )
  txt= paste(txt, tmp,fuki,goku, sep="" )
}

# 球面性
kazu=length(yui); sumi=0
mauW=MsW[c(NA,1,2,1,2,3,3)]; mauP=MsP[c(NA,1,2,1,2,3,3)]
if( kazu>0 )
{
  for(x in 1:kazu)
  {
    L= yui[x]
    if( L>1 )
    {
      shus= 0
      if( DF[3*(L%%2+1)]>1 | (L>5 & DF[3]+DF[6]>2) )
      {
        sumi=sumi+1; if(sumi==1) txt=paste(txt, "\n　有意性を示した自由度2以上の効果についてMauchlyの球面性検定を行った結果 (Table(tx9)参照)，", sep="" )
        if( is.nan(mauP[L])==1 ){ tmp="計算不能であった"; pch=NA; shus=shus+1 } else
        {
          pch= mauP[L]
          if( pch>1 ){ tmp="計算不能であった"; pch=NA; shus=shus+1 } else
            if( pch>0.05 ) tmp="有意でなかったことを確認した" else { tmp="有意であった"; shus=shus+1 }
        }
        if(sumi>1 & shus<1) goku="も" else goku="は"
        fuki= paste(" (Mauchly's W=",round(mauW[L],3),", p=",PV(pch),")。", sep="" )
        txt= paste(txt, NAM[L],"について",goku,tmp,fuki, sep="" )
        if( (pch<0.05 | is.na(pch)>0) & shus>0 )
        {
          goku="近似的に"
          if( sum(n==mean(n))>0 ) goku="" # unbalance
          if( shus==1 ) tmp=paste("このため",goku,"Greenhouse-Geisserの自由度調整係数 (ε) による修正検定を行った (Table(crP)参照)。結果として，", sep="")　else tmp="修正検定の結果，"
          txt= paste(txt,tmp, sep="" )
          if( L==3 | L==4 ) gyo=L-(L-3)*2 else gyo=L-1
          crp=crP[gyo,4]; fuki=""; ika=""
          if( crp<0.10 ){ tmp="であることを確認した"; if(crp>0.05) fuki="傾向" } else
          { tmp="性を得られなかった"; ika="以下，参考までに分析を進めることにする。" }
          txt= paste(txt, NAM[L],"は有意",fuki,tmp," (G-G corrected p=",PV(crp),")。",ika, sep="" )
        }
      }  } # if(L..)
  } # for(x)
} # if(kazu)
txt=paste(txt,"\n",sep="")

# 主効果
Pari=0; Mari=0; Tari=0
for( L in 1:3 ) #####
{
  pch= Pch[JUNB(L)]
  if( L==1 ){ levL=levA; hk=hkA[2,]; hyo=mcA; nn=hkA[1,] }
  if( L==2 ){ levL=levB; hk=hkB[2,]; hyo=mcB }
  if( L==3 ){ levL=levC; hk=hkC[2,]; hyo=mcC }
  nam= nama[L]

  if( pch<0.10 & min(Pch[Xari[L,]])>0.05 )
  {
    Pari= Pari+1
    if( Pari==1 ){ mata="\n　有意性を示した"; goku="" } else { mata=""; goku="は"; if(Pari==2) mata="また，" }
    txt= paste(txt, mata,NAM[L],"について",goku,"，", sep="" )
    hk=round(hk, 3 )
    if( levL==2 )
    {
      if( hk[1]>hk[2] ) tmp="大きい" else tmp="小さい"
      if( pch>0.05 ) goku="傾向がある" else goku=""
      txt= paste(txt, nam,"1の平均",hk[1],"が",nam,"2の平均",hk[2],"よりも有意に",tmp,goku,"ことが見いだされた。", sep="" )
    } else {
      Tari=1
      if(L<2) goku="プールドSDを用いた" else goku="対応のある"
      txt= paste(txt, goku,"ｔ検定による多重比較 (α=0.05, 両側検定) を行った結果，", sep="" )
      kazu= sum(hyo<0.10,na.rm=1 )
      if( kazu<1 ) txt=paste(txt,"どの２水準の平均の差も有意でなかった (adjusted ps>",PV(min(hyo,na.rm=1)),")。", sep="" ) else
      {
        ##
        sumi= 0
        for(x in 1:(levL-1)){ for(y in (x+1):levL){
          adP= hyo[ y-1, x ]
          if( adP<0.10 )
          {
            if(L==1)
            {
              tch= abs(hk[x]-hk[y])/( sqrt(1/nn[x]+1/nn[y])*plSD[1] )
              DFe= plSD[2]
            } else {
              dt= subset(AsBC,AsBC[,L+1]==x | AsBC[,L+1]==y )
              tx= t.test(dt[,5]~dt[,L+1],pair=1 )
              tch=abs(tx$stat); DFe=tx$para
            }
            keko="傾向がある"; if( adP<0.05 ) keko=""
            sumi= sumi+1
            if( sumi==2 ) mata="また" else mata=""
            if( sumi<kazu ){ tmp="こと"; goku="，" } else { tmp="ことが見いだされた"; goku="。" }
            fuki= paste(" t(",DFe,")=",round(tch,3)," adjusted p=",PV(adP),")",goku, sep="" )
            if( hk[x]>hk[y] ) mds="大きい" else mds="小さい"
            txt= paste(txt, mata,nam,x,"の平均",hk[x],"が",nam,y,"の平均",hk[y],"よりも有意に",mds,keko,tmp,fuki, sep="" )
          } # if(adP)
        }  } # for(x)(y)
        ##
      }
    } # if(levL)
  } # if(pch)
} #####

if( taju=="BH" ) TAJ="Benjamini & Hochberg (1995) "
if( taju=="holm" ) TAJ="Holm"
if( taju=="b" ) TAJ="Bonferroni"
if( Tari==1 ) txt=paste(txt, "\n　以上のp値の調整には",TAJ,"の方法を用いた。\n", sep="" ) else txt=paste(txt,"\n",sep="")
if( Tari==1 & taju=="BH" ) txt= paste(txt, "\n［引用文献］\n",
                                      "Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 58, 289-300.\n", sep="" )

txt= paste(txt,"\n",sep="" )
pch= min(Pch[c(4,7,9)] )
if( pch<0.10 ){ tmp="するか判断"; if( pch<0.05 ) tmp=""
if( Pch[10]>0.05 ) txt=paste(txt,"\n⇒js-STARの第二枠『一次の交互作用が有意のときに…』を実行",tmp,"してください。\n",sep="" ) }
if( Pch[10]<0.10 ){ tmp="するか判断"; if( Pch[10]<0.05 ) tmp=""
txt=paste(txt,"\n⇒js-STARの第三枠『二次の交互作用ＡxＢxＣが有意のときに…』を実行",tmp,"してください。\n",sep="" ) }
H=matrix(0,nr=1)
rownames(H)=""
colnames(H)=paste(" (多数回検定時のｐ値調整は",taju,"法による) ", sep="" )
nams= matrix(c("A","A","B","B","C","C"),nr=3)
TesA=function(水準=0,相手=0)
{
  tch= abs(hkA[2,水準]-hkA[2,相手])/( sqrt(1/hkA[1,水準]+1/hkA[1,相手])*plSD[1] )
  df2= plSD[2]
  cat(paste("t(",df2,")=",round(tch,4),", non-adjusted p=",round(pt(tch,df2,low=0)*2,4),"(両側)\n", sep="") )
}
TesB= function(水準=0,相手=0)
{
  kx= t.test(data~B,subset(AsBC,B==水準 | B==相手),pair=1)
  cat(paste("t(",kx$par,")=",round(kx$sta,4),", non-adjusted p=",round(kx$p.v,4),"(両側)\n", sep="") )
}
TesC=function(水準=1,相手=2)
{
  kx= t.test(data~C,subset(AsBC,C==水準 | C==相手),pair=1)
  cat(paste("t(",kx$par,")=",round(kx$sta,4),", non-adjusted p=",round(kx$p.v,4),"(両側)\n", sep="") )
}

options(digits=5)
options(scipen=5)
##############################
#  分散分析 ＡｓＢＣ-design
#  ３要因 混合計画Ⅱ
#
#  参加者間：要因Ａ
#  参加者内：要因Ｂ・Ｃ
H#############################
tx0 # 基本統計量（SD=不偏分散の平方根）
# 歪度,尖度,SEなどは『オプション』参照

tx1 # 分散分析ＡｓＢＣ
# 効果量ηp2 は偏イータ２乗

tx2 # 効果量ｆと検出力(1-β)
# 検出力0：水準間相関=0 (正負混在の場合)
# 検出力r：水準間相関=r (標本値から算出)

tx8 # 分散の均一性の検定（Bartlett Test）

tx9 # 球面性検定（df=1は不要）と自由度調整係数ε
crP # 球面性検定が有意のときの修正ｐ値
# ■参加者数Ｎ≦10ならH-F_ｐ参照可

hkA;mcA # 主効果Ａの平均と多重比較の調整後ｐ値

hkB;mcB # 主効果Ｂの平均と多重比較の調整後ｐ値

hkC;mcC # 主効果Ｃの平均と多重比較の調整後ｐ値

cat(txt) # 結果の書き方

# ■オプション：［↑］⇒行頭の♯を消す⇒［Enter］
# txD  # Median,歪度,尖度,SEなど
# ty2  # TypeⅡ_SS
# round(Bmat,3) # 要因Ｂの相関行列
# round(Cmat,3) # 要因Ｃの相関行列
# round(Tmat,3) # 要因Ｂ・Ｃの相関行列
# plSD # 主効果ＡのプールドSD
# TesA(水準=1, 相手=2 ) # 主効果Aの２水準ｔ検定
# TesB(水準=1, 相手=2 ) # 主効果Bの２水準ｔ検定
# TesC(水準=1, 相手=2 ) # 主効果Cの２水準ｔ検定
# write(txt,file.choose(),ap=T) # 結果のファイル保存
