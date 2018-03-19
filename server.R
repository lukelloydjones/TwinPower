require(shiny)
require(fpow)

shinyServer(function(input, output){ 
	
	values <- reactiveValues(default = 0)
    
    observeEvent(input$goButton,{
           values$default <- input$goButton
      })

	
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
    xncp_ml = ncparamF(2 * alpha, beta, 1, 10000000) # df2= 10000000 is an arbitrary 
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
    n_a_opt_ml = floor(xncp_ml / x_opt) # Optimum number of pairs
    n_a_mz_opt = round(n_a_opt_ml * p_a_ml)
    n_a_dz_opt = n_a_opt_ml - n_a_mz_opt
    # --------------------------------------
    incProgress(1/8, detail = paste("Calculated optimum number of pairs to detect A "))
    Sys.sleep(0.4)
    # --------------------------------------
    f_eb_mz = h2 + c2 + 1
    f_ew_mz = 1 - h2 - c2
    f_eb_dz = (h2 / 2) + c2 + 1 
    f_ew_dz = 1 - 0.5 * h2 - c2
    s11 <- seq(from = 0.01, to = 1, by = 0.01)
    s21 <- seq(from = 0.01, to = 1, by = 0.01)
    s1  <- expand.grid(s11, s21)   
    s1$fae <- c()
    b_mz = s1[, 1] + 2 * s1[, 2]
    w_mz = s1[, 1]
    b_dz = s1[, 1] + (3 / 2) * s1[, 2]  #not sure about this
    w_dz = s1[, 1] + (1 / 2) * s1[, 2]
    s1$fae = pmz * log(b_mz) + 
            pmz * log(w_mz) + 
            (1 - pmz) * log(b_dz) + 
            (1 - pmz) * log(w_dz) + 
            pmz * f_eb_mz / b_mz  +
            pmz * f_ew_mz / w_mz  + 
            (1 - pmz) * f_eb_dz / b_dz + 
            (1 - pmz) * f_ew_dz / w_dz 
    s1.min  = s1[s1$fae == min(s1$fae), ]
    x_ace   = pmz * (log(f_eb_mz) + log(f_ew_mz)) + 
              (1 - pmz) * (log(f_eb_dz) + log(f_ew_dz)) + 2
    xncp_ae = s1.min$fae - x_ace
    n_ae    = 1 + round((xncp_ml / xncp_ae))
    results.input <- c(pmz, s1.min$Var2, s1.min$Var1, s1.min$fae, xncp_ae, n_ae)    
    xncp_C  = (ndz + nmz) * xncp_ae
    power_C = pchisq(thres, 1, ncp = xncp_C, lower.tail = FALSE) 
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
      s$fae = p * log(b_mz) + 
              p * log(w_mz) + 
              (1 - p) * log(b_dz) + 
              (1 - p) * log(w_dz) + 
              p * f_eb_mz / b_mz  +
              p * f_ew_mz / w_mz  + 
              (1 - p) * f_eb_dz / b_dz + 
              (1 - p) * f_ew_dz / w_dz 
      s       = s[s$fae == min(s$fae),]
      x_ace   = p * (log(f_eb_mz) + log(f_ew_mz)) + 
               (1 - p) * (log(f_eb_dz) + log(f_ew_dz)) + 2
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
    # -------------------
    # First results table
    # ------------------- 
    output$view1 <- renderTable({
      gts  <- c("Central chi-squared threshold 1/2 Chi(1)", 
                "Required NCP for ML test (1-sided)",
                "Proportion of MZ twin pairs as defined by user")
      vals <- c(round(thres, 3),
                round(xncp_ml, 3),
                round(pmz, 3))
      df.1  <- data.frame(gts, vals)
      colnames(df.1) <- c(" ", " ")
      df.1
    }, digits = 3, include.rownames = FALSE)
    # -------------------
    # Second results table
    # -------------------
    output$view2 <- renderTable({
      gts  <- c("Power to detect A given input parameters",
                "Number of pairs to detect A for given ratio of MZ pairs",
                "Optimum proportion of MZ pairs",
                "Optimum number of pairs to achieve required power")
      vals <- c(round(powerA, 3),
                format(n_a_ml, scientific = FALSE),
                round(p_a_ml, 3),
                format(n_a_opt_ml, scientific = FALSE))
      df.1  <- data.frame(gts, vals)
      colnames(df.1) <- c(" ", " ")
      df.1
    }, digits = 3, include.rownames = FALSE)
    # -------------------
    # Third results table
    # -------------------
    output$view3 <- renderTable({
    gts  <- c("Power to detect C given input parameters",
              "Number of pairs to detect C for given MZ ratio",
              "Optimum MZ ratio to achieve required power",
              "Minimum sample size to achieve required power")
    vals <- c(round(power_C, 3),
              format(results.input[6], scientific = FALSE),
              out1[1],
              format(out1[6], scientific = FALSE))
      df.1  <- data.frame(gts, vals)
      colnames(df.1) <- c(" ", " ")
      df.1
      }, digits = 3, include.rownames = FALSE)
      HTML("")
  })
  })
  
  getVals2 <- eventReactive(1, {
    # -------------------
    # First results table
    # ------------------- 
    output$view1 <- renderTable({
      gts  <- c("Central chi-squared threshold 1/2 Chi(1)", 
                "Required NCP for ML test (1-sided)",
                "Proportion of MZ twin pairs as defined by user")
      vals <- c(2.706,
                6.182,
                0.50)
      df.1  <- data.frame(gts, vals)
      colnames(df.1) <- c(" ", " ")
      df.1
    }, digits = 3, include.rownames = FALSE)
    # -------------------
    # Second results table
    # -------------------
    output$view2 <- renderTable({
      gts  <- c("Power to detect A given input parameters",
                "Number of pairs to detect A for given ratio of MZ pairs",
                "Optimum proportion of MZ pairs",
                "Optimum number of pairs to achieve required power")
      vals <- c(0.659,
                format(60, scientific = FALSE),
                0.572,
                as.integer(57))
      df.1  <- data.frame(gts, vals)
      colnames(df.1) <- c(" ", " ")
      df.1
    }, digits = 3, include.rownames = FALSE)
    # -------------------
    # Third results table
    # -------------------
    output$view3 <- renderTable({
    gts  <- c("Power to detect C given input parameters",
              "Number of pairs to detect C for given MZ ratio",
              "Optimum MZ ratio to achieve required power",
              "Minimum sample size to achieve required power")
    vals <- c(0.157,
              format(725, scientific = FALSE),
              0.140,
              as.integer(488))
      df.1  <- data.frame(gts, vals)
      colnames(df.1) <- c(" ", " ")
      df.1
      }, digits = 3, include.rownames = FALSE)
      HTML("")
  })
  
  output$text1 <- renderUI({
   	 if(values$default == 0)
     {
     	getVals2()
     } else {
  	    getVals()
  	 }
   })
   
 
})

    
