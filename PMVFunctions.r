# Compute PMV with different algorithms -----------------------------------

# Functions to compute hygrothermal values
# Florian Antretter
# Date: 2018-04-12
# Version: 0.1


# Algorithm used in old R -------------------------------------------------

calcPMV_IBP <- function(ta = 20, tr = 20, rh = 50, va = 0.1, met = 1.2, clo = 0.7, W = 0) {
  ### Predicted Mean Vote Calculation according to ISO 7730
  ### from  
  ### ta = ambient temperature in °C
  ### top = operative temperature in °C
  ### rh = relative humidity in %
  ### va = air velocity in m/s
  ### met = metabolic rate
  ### clo = clothing insulation value
  ### W rate of mechanical work accomplishes (usually 0)
  
  PMV <- c()
  # compute metabolism and clothing thermal resistance
  M     <- met * 58             # metabolic rate [W/m2]
  Icl   <- clo * 0.155          # thermal resistance of clothing [m2K/W]
  
  # compute partial water vapour pressure [Pa] according to ISO 7726
  pa    <- c() 
  pa    <- rh*(6.11*exp(17.27*ta/(237.3+ta)))  
  
  # compute clothing area coefficient [-]
  fcl   <- c()   
  fcl   <- ifelse(Icl <= 0.078, 1+1.29*Icl, 1.05+0.645*Icl)
  
  # compute mean radiant temperature [°C] from operative temperature [°C]
  #tr   <- c()                                                   # mean radiant temperature [?C]
  #tr   <- top * 2 - ta
  
  # compute clothing surface temperature [°C]
  tcl <- c()                                                 # surface temperature of clothing [?C]
  for(i in 1:length(ta)) {
    tcl.function <- function(tcl.i) {
      (tcl.i - (35.7 - 0.028*(M[i]-W[i]) - Icl[i]*(3.96e-8*fcl[i]*((tcl.i+273)^4-(tr[i]+273)^4)+fcl[i]*max(2.38*(abs(tcl.i-ta[i]))^0.25,12.1*sqrt(va[i]))*(tcl.i-ta[i]))))
    }
    
    tcl.i  <- uniroot(tcl.function,lower = -50, upper = 70)
    tcl[i] <- tcl.i$root
  }
  
  # compute convective heat transfer coefficient [-]
  hc   <- c()
  hc   <- max(2.38*(abs(tcl-ta))^0.25,12.1*sqrt(va))
  #for (i in 1:length(ta))                                 
  #{hc[i] <- max(2.38*(abs(tcl[i]-ta[i]))^0.25,12.1*sqrt(va[i]))}
  
  ### output: predicted mean vote [-]
  (0.303*exp(-0.036*M)+0.028) * ((M-W)-0.00305*(5733-6.99*(M-W)-pa)-0.42*((M-W)-58.15)-1.7e-5*M*(5867-pa)-0.0014*M*(34-ta)
                                        -3.96e-8*fcl*((tcl+273)^4-(tr+273)^4)-fcl*hc*(tcl-ta))
  
  #return(PMV)
}


# Calc PMV according to the code implemented in WP ------------------------

require(Rcpp)
cppFunction(code='
            double calcPMV_WP(double t_a, double t_r, double rh, double vel, double met, double clo, double wme)
            {
            double t_aa, t_ra,
            ps, pa, I_cl, M, W, H, f_cl, h_cf, h_cn, h_c,
            eps, p1, p2, p3, p4, p5, xn, xf,
            n, n_max,
            t_cla, t_cl,
            hl1, hl2, hl3, hl4, hl5, hl6,
            ts, result;
            
            //initialize
            eps = 0.00015;
            n_max = 150;
            
            //saturated vapour pressure [kPa]
            ps = exp(16.6536 - 4030.183 / (t_a + 235));
            
            //water vapor pressure [Pa]
            pa = rh * 1000 * ps;
            
            //thermal insulation of the clothing [m²K/W]
            I_cl = 0.155 * clo;
            
            //metabolic rate [W/m²]
            M = met * 58.15;
            
            //external work [W/m²]
            W = wme * 58.15;
            
            //internal heat production in the human body
            H = M - W;
            
            //clothing area factor
            if (I_cl < 0.078) f_cl = 1 + 1.29 * I_cl;
            else   f_cl = 1.05 + 0.645 * I_cl;
            //heat transfer coefficient by forced convection
            h_cf = 12.1 * sqrt(vel);
            //air temperature [K]
            t_aa = t_a + 273;
            //mean radiant temperature [K]
            t_ra = t_r + 273;
            //CALCULATE SURFACE TEMPERATURE OF CLOTHING BY ITERATION
            t_cla = t_aa + (35.5 - t_a) / (3.5 * (6.45 * I_cl + 0.1));
            //initial guess for surface temperature of clothing
            n = 0;
            p1 = I_cl * f_cl;
            p2 = p1 * 3.96;
            p3 = p1 * 100;
            p4 = p1 * t_aa;
            p5 = 308.7 - 0.028 * H + p2 * pow(0.01 * t_ra, 4);  //(0.01 * t_ra) ^ 4;
            xn = 0.01 * t_cla;
            xf = xn;
            do  //Do
            {
            n = n + 1;
            //If n > 150 Then  Exit Do  End If
            if (n < 150) 
            {
            xf = 0.5 * (xf + xn);
            //heat transfer coefficient by natural convection
            h_cn = 2.38 * pow(abs(100 * xf - t_aa),0.25);  // h_cn = 2.38 * Abs(100 * xf - t_aa) ^ 0.25
            
            if (h_cf > h_cn) h_c = h_cf;
            else  h_c = h_cn;
            
            xn = (p5 + p4 * h_c - p2 * pow(xf,4)) / (100 + p3 * h_c);
            }
            }
            //Loop While (Abs(xn - xf) > eps)
            while ((abs(xn - xf) > eps) && (n <= 150));
            
            if (n < 150) 
            {
            //surface temperature of clothing
            t_cl = 100 * xn - 273;
            
            //HEAT LOSS COMPONENTS
            //heat loss diff. through skin
            hl1 = 3.05 * 0.001 * (5733 - 6.99 * H - pa);
            
            //heat loss by sweating (comfort)
            if (H > 58.15) hl2 = 0.42 * (H - 58.15);
            else  hl2 = 0;
            //End If
            
            //latent respiration heat loss
            hl3 = 1.7 * 0.00001 * M * (5867 - pa);
            
            //dry respiration heat loss
            hl4 = 0.0014 * M * (34 - t_a);
            
            //heat loss by radiation
            hl5 = 3.96 * f_cl * (pow(xn, 4) - pow((t_ra * 0.01), 4));
            
            //heat loss by convection
            hl6 = f_cl * h_c * (t_cl - t_a);
            
            //CALCULATE PMV
            //thermal sensation trans. coefficient
            ts = 0.303 * exp(-0.036 * M) + 0.028;
            //predicted mean vote
            //PMV_ISO = ts * (H - hl1 - hl2 - hl3 - hl4 - hl5 - hl6);
            result = ts * (H - hl1 - hl2 - hl3 - hl4 - hl5 - hl6);
            }
            else
            //PMV_ISO = 9999;
            result = 0; // ??  9999
            return result;
            }
            ')


# Calc PMV according to an adjusted version of the implementation  --------

calcPMV <- function(ta, tr, rh, vel, met=1, clo=.5, wme=0, basMet=58.15){
  
  m   <- met * basMet 
  w   <- wme * basMet
  mw  <- m - w
  icl <- .155 * clo
  pa  <- rh * 10 * exp(16.6536 - (4030.183 / (ta + 235)))
  
  # Compute the corresponding fcl value
  fcl <- ifelse(icl <= .078, 1 + 1.29 * icl, 1.05 + .645 * icl)
  
  fcic <- icl * fcl
  p2   <- fcic * 3.96
  p3   <- fcic * 100
  tra  <- tr + 273
  taa  <- ta + 273
  p1   <- fcic * taa
  p4   <- 308.7 - .028 * mw + p2 * (tra / 100) ^ 4
  
  # First guess for surface temperature
  tclA <- taa + (35.5-ta) / (3.5 * (6.45 * icl + .1))
  xn   <- tclA / 100
  xf   <- xn
  hcf  <- 12.1 * (vel) ^ .5
  noi  <- 0
  eps  <- .00015
  
  
  # COmPUTE sURFAce TEmPEratuRE OF cloTHING BY ITEraTIONs
  while (noi < 200){
    xf  <- (xf + xn) / 2
    hcn <- 2.38 * abs(100 * xf - taa) ^ .25
    hc  <- ifelse (hcf > hcn, hcf, hcn)
    xn  <- (p4 + p1 * hc - p2 * xf ^ 4) / (100 + p3 * hc)
    noi <- noi + 1
    #if(noi > 1 & abs(xn - xf) <= eps){break} 
  }
  tcl <- 100 * xn - 273
  
  # COmPUTE pmv
  
  pm1 <- 3.96 * fcl * (xn ^ 4 - (tra / 100) ^ 4)
  pm2 <- fcl * hc * (tcl - ta)
  pm3 <- .303 * exp(-.036 * m) + .028
  pm4 <- ifelse (mw > basMet, .42 * (mw - basMet), 0)
  pm5 <- 3.05 * .001 * (5733 - 6.99 * mw - pa)
  pm6 <- 1.7 * .00001 * m * (5867 - pa) + .0014 * m * (34 - ta)
  pm3 * (mw - pm5 - pm4 - pm6 - pm1 - pm2)
  
  #ppd <- 100 - 95 * exp(-.03353 * pmv ^ 4 - .2179 * pmv ^ 2)
  #data.frame(pmv, ppd)
  
  #return(pmv)
}
