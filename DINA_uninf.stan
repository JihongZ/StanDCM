    
data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
}
   
parameters{
  simplex[Nc] Vc;
   real<lower=0> l1_12 ;
   real<lower=0> l2_223 ;
   real<lower=0> l3_213 ;
   real<lower=0> l4_213 ;
   real<lower=0> l5_212 ;
   real<lower=0> l6_212 ;
   real<lower=0> l7_212 ;
   real<lower=0> l8_223 ;
   real<lower=0> l9_213 ;
   real l1_0 ;
   real l2_0 ;
   real l3_0 ;
   real l4_0 ;
   real l5_0 ;
   real l6_0 ;
   real l7_0 ;
   real l8_0 ;
   real l9_0 ;
 }
 
 transformed parameters{
 matrix[Ni, Nc] PImat;
 vector[Ni] gParm;
 vector[Ni] sParm;
  real l2_12 ;
   real l2_13 ;
   real l3_11 ;
   real l3_13 ;
   real l4_11 ;
   real l4_13 ;
   real l5_11 ;
   real l5_12 ;
   real l6_11 ;
   real l6_12 ;
   real l7_11 ;
   real l7_12 ;
   real l8_12 ;
   real l8_13 ;
   real l9_11 ;
   real l9_13 ;
   l2_12 =0 ;
   l2_13 =0 ;
   l3_11 =0 ;
   l3_13 =0 ;
   l4_11 =0 ;
   l4_13 =0 ;
   l5_11 =0 ;
   l5_12 =0 ;
   l6_11 =0 ;
   l6_12 =0 ;
   l7_11 =0 ;
   l7_12 =0 ;
   l8_12 =0 ;
   l8_13 =0 ;
   l9_11 =0 ;
   l9_13 =0 ;
   gParm[1]=inv_logit(l1_0);
  gParm[2]=inv_logit(l2_0);
  gParm[3]=inv_logit(l3_0);
  gParm[4]=inv_logit(l4_0);
  gParm[5]=inv_logit(l5_0);
  gParm[6]=inv_logit(l6_0);
  gParm[7]=inv_logit(l7_0);
  gParm[8]=inv_logit(l8_0);
  gParm[9]=inv_logit(l9_0);
  sParm[1]=1-inv_logit(l1_0+l1_12);
  sParm[2]=1-inv_logit(l2_0+l2_12+l2_13+l2_223);
  sParm[3]=1-inv_logit(l3_0+l3_11+l3_13+l3_213);
  sParm[4]=1-inv_logit(l4_0+l4_11+l4_13+l4_213);
  sParm[5]=1-inv_logit(l5_0+l5_11+l5_12+l5_212);
  sParm[6]=1-inv_logit(l6_0+l6_11+l6_12+l6_212);
  sParm[7]=1-inv_logit(l7_0+l7_11+l7_12+l7_212);
  sParm[8]=1-inv_logit(l8_0+l8_12+l8_13+l8_223);
  sParm[9]=1-inv_logit(l9_0+l9_11+l9_13+l9_213);
  PImat[1,1]=inv_logit(l1_0);
  PImat[2,1]=inv_logit(l2_0);
  PImat[3,1]=inv_logit(l3_0);
  PImat[4,1]=inv_logit(l4_0);
  PImat[5,1]=inv_logit(l5_0);
  PImat[6,1]=inv_logit(l6_0);
  PImat[7,1]=inv_logit(l7_0);
  PImat[8,1]=inv_logit(l8_0);
  PImat[9,1]=inv_logit(l9_0);
  PImat[1,2]=inv_logit(l1_0);
  PImat[2,2]=inv_logit(l2_0+l2_13);
  PImat[3,2]=inv_logit(l3_0+l3_13);
  PImat[4,2]=inv_logit(l4_0+l4_13);
  PImat[5,2]=inv_logit(l5_0);
  PImat[6,2]=inv_logit(l6_0);
  PImat[7,2]=inv_logit(l7_0);
  PImat[8,2]=inv_logit(l8_0+l8_13);
  PImat[9,2]=inv_logit(l9_0+l9_13);
  PImat[1,3]=inv_logit(l1_0+l1_12);
  PImat[2,3]=inv_logit(l2_0+l2_12);
  PImat[3,3]=inv_logit(l3_0);
  PImat[4,3]=inv_logit(l4_0);
  PImat[5,3]=inv_logit(l5_0+l5_12);
  PImat[6,3]=inv_logit(l6_0+l6_12);
  PImat[7,3]=inv_logit(l7_0+l7_12);
  PImat[8,3]=inv_logit(l8_0+l8_12);
  PImat[9,3]=inv_logit(l9_0);
  PImat[1,4]=inv_logit(l1_0+l1_12);
  PImat[2,4]=inv_logit(l2_0+l2_12+l2_13+l2_223);
  PImat[3,4]=inv_logit(l3_0+l3_13);
  PImat[4,4]=inv_logit(l4_0+l4_13);
  PImat[5,4]=inv_logit(l5_0+l5_12);
  PImat[6,4]=inv_logit(l6_0+l6_12);
  PImat[7,4]=inv_logit(l7_0+l7_12);
  PImat[8,4]=inv_logit(l8_0+l8_12+l8_13+l8_223);
  PImat[9,4]=inv_logit(l9_0+l9_13);
  PImat[1,5]=inv_logit(l1_0);
  PImat[2,5]=inv_logit(l2_0);
  PImat[3,5]=inv_logit(l3_0+l3_11);
  PImat[4,5]=inv_logit(l4_0+l4_11);
  PImat[5,5]=inv_logit(l5_0+l5_11);
  PImat[6,5]=inv_logit(l6_0+l6_11);
  PImat[7,5]=inv_logit(l7_0+l7_11);
  PImat[8,5]=inv_logit(l8_0);
  PImat[9,5]=inv_logit(l9_0+l9_11);
  PImat[1,6]=inv_logit(l1_0);
  PImat[2,6]=inv_logit(l2_0+l2_13);
  PImat[3,6]=inv_logit(l3_0+l3_11+l3_13+l3_213);
  PImat[4,6]=inv_logit(l4_0+l4_11+l4_13+l4_213);
  PImat[5,6]=inv_logit(l5_0+l5_11);
  PImat[6,6]=inv_logit(l6_0+l6_11);
  PImat[7,6]=inv_logit(l7_0+l7_11);
  PImat[8,6]=inv_logit(l8_0+l8_13);
  PImat[9,6]=inv_logit(l9_0+l9_11+l9_13+l9_213);
  PImat[1,7]=inv_logit(l1_0+l1_12);
  PImat[2,7]=inv_logit(l2_0+l2_12);
  PImat[3,7]=inv_logit(l3_0+l3_11);
  PImat[4,7]=inv_logit(l4_0+l4_11);
  PImat[5,7]=inv_logit(l5_0+l5_11+l5_12+l5_212);
  PImat[6,7]=inv_logit(l6_0+l6_11+l6_12+l6_212);
  PImat[7,7]=inv_logit(l7_0+l7_11+l7_12+l7_212);
  PImat[8,7]=inv_logit(l8_0+l8_12);
  PImat[9,7]=inv_logit(l9_0+l9_11);
  PImat[1,8]=inv_logit(l1_0+l1_12);
  PImat[2,8]=inv_logit(l2_0+l2_12+l2_13+l2_223);
  PImat[3,8]=inv_logit(l3_0+l3_11+l3_13+l3_213);
  PImat[4,8]=inv_logit(l4_0+l4_11+l4_13+l4_213);
  PImat[5,8]=inv_logit(l5_0+l5_11+l5_12+l5_212);
  PImat[6,8]=inv_logit(l6_0+l6_11+l6_12+l6_212);
  PImat[7,8]=inv_logit(l7_0+l7_11+l7_12+l7_212);
  PImat[8,8]=inv_logit(l8_0+l8_12+l8_13+l8_223);
  PImat[9,8]=inv_logit(l9_0+l9_11+l9_13+l9_213);
}
 
model {
    vector[Nc] contributionsC;
    vector[Ni] contributionsI;

    //Prior
    l1_12~normal(0,15);
                    l2_223~normal(0,15);
    l3_213~normal(0,15);
    l4_213~normal(0,15);
    l5_212~normal(0,15);
    l6_212~normal(0,15);
    l7_212~normal(0,15);
    l8_223~normal(0,15);
    l9_213~normal(0,15);
    l1_0~normal(0,15);
    l2_0~normal(0,15);
    l3_0~normal(0,15);
    l4_0~normal(0,15);
    l5_0~normal(0,15);
    l6_0~normal(0,15);
    l7_0~normal(0,15);
    l8_0~normal(0,15);
    l9_0~normal(0,15);
    Vc~dirichlet(rep_vector(2.0, Nc)); 
  

  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
        else
          contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
      }
      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
    }
  target+=log_sum_exp(contributionsC);
  }
   
}  
  

generated quantities {
  vector[Ni] log_lik[Np];
  vector[Ni] contributionsI;
  matrix[Ni,Nc] contributionsIC;
  //Posterior
  for (iterp in 1:Np){
    for (iteri in 1:Ni){
      for (iterc in 1:Nc){
        if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
        else
          contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
        contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
      }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
  }
}
  