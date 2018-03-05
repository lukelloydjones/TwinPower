require(shiny)
require(fpow)

shinyServer(function(input, output){
  

    
    getVals <- eventReactive(input$goButton, {

    nmz = input$nmz
    ndz = input$ndz
    Va  = input$Va
    Vc  = input$Vc
    Ve  = input$Ve
    alpha = input$alpha
    power = input$power
    # ------------------------------------------------------------
    # Begin prograss bar
    # ------------------------------------------------------------
    withProgress(message = 'Performing calculations', value = 0, {
    # Initial calculations
    Vp = Va + Vc + Ve
    h2 = Va/Vp
    c2 = Vc/Vp
    t_mz = h2 + c2
    t_dz = (0.5 * h2) + c2
    pmz = nmz / (nmz + ndz)
    t_ave = pmz * t_mz + (1-pmz) * t_dz
    beta = 1 - power
    # --------------------------------------
    incProgress(1/8, detail = paste("Initial calculations complete"))
    Sys.sleep(0.4)
    # --------------------------------------
    # Calculate threshold (null) and required NCP, assuing a 50:50 0:chi(1) 
    # distribution under the null
    # --------------------------------------
    thres   = qchisq(1 - (2 * alpha), 1)
    xncp_ml = ncparamF(2 * alpha, beta, 1, 10000) # df2= 10000 is an arbitrary 
                                             # large number so that F(1,1000) ~ chi(1)
    # --------------------------------------
    incProgress(1/8, detail = paste("Calculated threshold and required NCP "))
    Sys.sleep(0.4)
    # --------------------------------------
    # Maximum Likelihood: ACE vs CE (i.e. detection of A)
    # --------------------------------------
    x_detA = log(1 - t_ave ^ 2) - pmz * log(1 - t_mz ^ 2) - (1-pmz) * log(1 - t_dz ^ 2)
    xncp_A = (ndz + nmz) * x_detA
    powerA = pchisq(thres, 1, ncp = xncp_A, lower.tail = FALSE) # The power to detect A
    n_a_ml = 1 + round(xncp_ml / x_detA )
    n_a_mz = round(n_a_ml * pmz)
    n_a_dz = n_a_ml - n_a_mz
    # --------------------------------------
    incProgress(1/8, detail = paste("Calculated ACE and CE maximum likelihood "))
    Sys.sleep(0.4)
    # --------------------------------------
    # Optimum number of pairs to detect A
    delta = t_mz - t_dz
    con   = log(1 - t_mz ** 2) - log(1 - t_dz ** 2)
    opt_t_ave = (delta - sqrt(delta ** 2 + con ** 2) ) / con
    p_a_ml    = (opt_t_ave - t_dz) / delta
    t_ave_opt =  p_a_ml * t_mz + (1-p_a_ml) * t_dz
    x_opt     = log(1 - t_ave_opt ** 2) - p_a_ml * log(1 - t_mz ** 2) - 
                   (1-p_a_ml) * log(1 - t_dz ** 2)
    n_a_opt_ml = round(xncp_ml / x_opt) # Optimum number of pairs
    n_a_mz_opt = round(n_a_opt_ml * p_a_ml)
    n_a_dz_opt = n_a_opt_ml - n_a_mz_opt
    # --------------------------------------
    incProgress(1/8, detail = paste("Calculated optimum number of pairs to detect A "))
    Sys.sleep(0.4)
    # --------------------------------------
    # MAXIMUM LIKELIHOOD: ACE vs AE (i.e. detection of C)
    # From Appendix A of Visscher 2004 (Twin Research)
    f_ew_mz = 1 - h2 - c2
    f_eb_mz = f_ew_mz + 2 * (h2 + c2)
    f_ew_dz = 1 - 0.5 * h2 - c2
    f_eb_dz = f_ew_dz + 2 * (0.5 * h2 + c2)
    vt_mz = (1 - t_mz ** 2) ** 2 / pmz
    vt_dz = (1 - t_dz ** 2) ** 2 / (1 - pmz)
    h2_ae = (t_mz / vt_mz + 2 * t_dz / (4 * vt_dz)) / (1 / vt_mz + 1 / (4 * vt_dz))
    x_ae  = (pmz * log(1 + h2_ae) + pmz * log(1 - h2_ae) + (1 - pmz) * log(1 + 0.5 * h2_ae) 
            + (1 - pmz) * log(1 - 0.5 * h2_ae) + 
            pmz * f_eb_mz / (1 + h2_ae) + pmz * f_ew_mz / (1 - h2_ae) +
            (1 - pmz) * f_eb_dz / (1 + 0.5 * h2_ae) + (1 - pmz) * f_ew_dz / (1 - 0.5 * h2_ae))
    x_ace = pmz * (log(f_eb_mz) + log(f_ew_mz)) + (1 - pmz) * (log(f_eb_dz) + log(f_ew_dz) ) + 2.0
    xncp_ae = x_ae - x_ace
    xncp_C = (ndz+nmz) * xncp_ae
    power_C = pchisq(thres, 1, ncp = xncp_C, lower.tail = FALSE) # The power to detect C
    n_ae = 1 + round(xncp_ml / xncp_ae)
    n_c_mz = round( pmz * n_ae )
    n_c_dz = n_ae - n_c_mz
    # --------------------------------------
    incProgress(1/8, detail = paste("Calculated optimum maximum likelihood for ACE vs AE model"))
    Sys.sleep(0.4)
    # --------------------------------------
    # --------------------------------------
    incProgress(1/8, detail = paste("Performing grid search"))
    Sys.sleep(0.4)
    # --------------------------------------
    # for different values of p grid search over different values of sige and siga
    out1 <- matrix(0, 50, 6)
    pp   <- seq(0.01, 0.5, 0.01)  
    for (j in 1:length(pp)){
      #print(j)  
      p  <- pp[j]
      s1 <- seq(from = 0.01, to = 1, by = 0.01)
      s2 <- seq(from = 0.01, to = 1, by = 0.01)
      s  <- expand.grid(s1, s2)   
      s$fae <- c()
      b_mz = s[, 1] + 2 * s[, 2]
      w_mz = s[, 1]
      b_dz = s[, 1] + (3 / 2) * s[, 2]  #not sure about this
      w_dz = s[, 1] + (1 / 2) * s[, 2]
      s$fae = p * log(b_mz) + p * log(w_mz) + (1 - p) * log(b_dz) + 
                  (1 - p) * log(w_dz) + p * f_eb_mz / b_mz + p * f_ew_mz / w_mz + 
                  (1-p) * f_eb_dz / b_dz + (1 - p) * f_ew_dz / w_dz 
      s       = s[s$fae == min(s$fae),]
      x_ace   = p * (log(f_eb_mz) + log(f_ew_mz)) + (1 - p) * (log(f_eb_dz) + log(f_ew_dz)) + 2
      xncp_ae = s$fae - x_ace
      n_ae    = 1 + round((xncp_ml / xncp_ae))
      out1[j,] <- (c(p, s$Var2, s$Var1, s$fae, xncp_ae, n_ae))
    }
    out1 <- matrix(out1[out1[, 6] == min(out1[, 6]), ], ncol = 6)
    if (dim(out1)[1] > 1)
    {
      out1 <- out1[out1[, 1] == min(out1[, 1]), ]
    }
    # --------------------------------------
    incProgress(1/8, detail = paste("Grid search complete "))
    Sys.sleep(0.4)
    # --------------------------------------
    # --------------------------------------
    incProgress(1/8, detail = paste("Outputting results"))
    Sys.sleep(0.4)
    # --------------------------------------
    #Output/Results
    HTML(paste(
        "RESULTS", '<br/>',
        "Central chisquare threshold 1/2 Chi(1)", round(thres, 4), '<br/>',
        "Required NCP for ML test (1-sided)", round(xncp_ml,4), '<br/>',
        "***** MAXIMUM LIKELIHOOD ACE vs. CE model ****",  '<br/>',
        "ML: NCP per pair for given pmz", round(x_detA,4), '<br/>',
        "Power to detect A given input parameters", round(powerA,3), '<br/>',
        "Number of pairs to detect A for given pmz", round(n_a_ml),  '<br/>',
        "ML: Optimum proportion of MZ", round(p_a_ml), '<br/>',
        "Optimal number of pairs", n_a_opt_ml,  '<br/>',
        "***** MAXIMUM LIKELIHOOD ACE vs. AE model ****",  '<br/>',
        "A", round(h2_ae, 4), '<br/>',
        "E", round(1-h2_ae, 4), '<br/>',
        "-2ML(AE)", round(x_ae, 4), '<br/>',
        "LRT", round(xncp_ae, 4), '<br/>',
        "Power to detect C given input parameters", round(power_C,3),  '<br/>',
        "Number of pairs to detect C for given pmz", round(n_ae),  '<br/>',
        "Best MZ ratio", out1[1], '<br/>',
        "Minimum sample size", out1[6])
       ) 
  })
  })
    output$text1 <- renderUI({
  	  getVals()
    })
})

    
