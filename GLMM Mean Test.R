# rm(list = ls()) #FOR A 'CLEAN' RUN, RESTART Rstudio: rm(list = ls()) is unfashionable
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2021 03 01
#     MODIFIED:	James Foster              DATE: 2021 07 22
#
#  DESCRIPTION: Loads a text files, fits a mixed effects logistic model and
#               estimates the p-value for mean ≤ reference mean.
#               
#       INPUTS: A ".csv" table with columns for experiment type ("Test"),
#               individual ID ("Animal_name"), 
#               proportion of correct choices ("%"),total number of choices ("total Choice number"), NOT IMPLEMENTED
#               correct choices ("correct_choices"), incorrect choices ("incorrect_choices").
#               User should specify h0 reference correct choice rate (line 50).
#               
#      OUTPUTS: Plots of confidence intervals (.pdf).
#
#	   CHANGES: - two-tailed test
#             - z-score used to calculate p-value
#             - user-defined column names
#             - two-tailed cumulative probability plotted as calculated
#
#   REFERENCES: Bates D, Maechler M, Bolker B, Walker S (2015).
#               Fitting Linear Mixed-Effects Models Using lme4.
#               Journal of Statistical Software, 67(1), 1-48.
#               doi:10.18637/jss.v067.i01.
#
#    EXAMPLES:  Fill out user input (lines 50-56), then press ctrl+shift+s to run
# 
#TODO   ---------------------------------------------
#TODO   
#- Read in data   +
#- Fit models     +
#- Extract fits   + 
#- Calculate p-value  +
#- Plot CI        +
#- Calculate z score +
#- Two-tailed test   +       
#- Worked examples
#- Back-calculate from % correct
#- Save PDFs

# Useful functions --------------------------------------------------------
# . Load package ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
require(lme4)#package for fitting GLMMs (generalised linear mixed-effects models)

# Input Variables ----------------------------------------------------------
#  .  User input -----------------------------------------------------------
h0 = 0.5#reference choice rate
h1 = 'two_tailed'#alternative hypothesis, mean "greater" than reference or "less", or either "two_tailed" [sic]
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
individual_factor = "Animal_name" #The title of the column; NO SPACES PLEASE
correct_name = "correct_choices" #The title of the column with positive choices; NO SPACES PLEASE
incorrect_name = "incorrect_choices" #The title of the column with negative choices; NO SPACES PLEASE

#Check the operating system and assign a logical flag (T or F)
sys_win <- Sys.info()[['sysname']] == 'Windows'
#On computers set up by JMU Würzburg, use user profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp <- gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp <- Sys.getenv('HOME')#Life was easier on Mac
}

# Simulate data (not used) ------------------------------------------------
  # n_animals = 10
  # n_choices = 20
  # correct_sim = rbinom(n = n_animals,
  #                      size = n_choices,
  #                      prob = 0.40
  #                      )
  # sim = data.frame(Animal_name = 1:n_animals,
  #                  correct_choices = correct_sim,
  #                  incorrect_choices = n_choices - correct_sim
  #                  )
  # write.table(x = sim,
  #             file = file.path(ltp,'Documents', "simulated_data.csv"),
  #             sep = csv_sep,
  #             row.names = FALSE
  #             )

# . Select files ---------------------------------------------------------

# set path to files
if(sys_win){#choose.files is only available on Windows
  message('\n\nPlease select the ".csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file  <- choose.files(
    default = file.path(ltp,'Documents', "*.csv"),#For some reason this is not possible in the "root" user
    caption = 'Please select the ".csv" file'
  )
}else{
  message('\n\nPlease select the ".csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file <- file.choose(new=F)
}
#show the user the path they have selected
if(is.null(path_file))
  {stop('No file selected.')}else
  {print(path_file)}


# Read in file ------------------------------------------------------------
mdata = read.table(file = path_file,#read from user-selected file
                    header = T,#read the file header to use for variable names
                    sep = csv_sep#,#values are separated by the user-specified character
                    #other parameters can be added here for troubleshooting
                    )

View(mdata)#show the user the data that was
#TODO, make this more flexible
mdata$ID<- mdata[,individual_factor]#column contains ID in some form
mdata$correct_choices<- mdata[,correct_name]#column contains total positive per animal
mdata$incorrect_choices<- mdata[,incorrect_name]#column contains total negative per animal
# Fit model ---------------------------------------------------------------
frm.1 = formula(#response variable is a pair of columns (correct, incorrect)
          cbind(correct_choices, incorrect_choices) ~ #responses are binomial choice counts
                1 + #the predictor is a single mean,
                (1|ID) # individuals have their own random-effects means
               )
mod.1 = glmer(formula = frm.1,#fit this formula to a _generalised_ linear model
               data = mdata,#using this data (i.e. header names match formula variables)
               family = binomial(link = 'logit')#mean and s.d. are fitted on logit scale https://en.wikipedia.org/wiki/Logit 
               )
#N.B. this is the simplest possible model including random effects,
#so no model selection is necessary

#inspect model
summary(mod.1)
#AIC:   Akaike Information Criterion, lower = better
#Random effects; Std. Dev.: should be >0 or individual effects were not fitted
#Fixed effects; Estimate: logit(mean success rate), 1/(1+exp(-Estimate)) = mean success rate

#N.B. for two tailed tests against a reference choice rate of 50%
#the p value from summary(mod.1) should be the same as the one printed
#to the command line. The null hypothesis for this test is that 
#log(correct / incorrect) = 0, which is the same as the probability being 50%
#since log(50/50) = log(1) = 0.

#  .  Derive variables	----------------------------------------------

#transform fixed effects mean to odds, logit is equivalent to log(odds)
message('correct choices per incorrect choice\n',
        round( exp(fixef(mod.1)), 3)
        )#X correct choices for every incorrect choice
#transform fixed effects mean to logistic probability (proportion correct)
mod.mu = fixef(mod.1)
message('mean correct choice rate\n',
        round(plogis(mod.mu), 3)*100,
        '%')
#extract fixed-effects standard deviation
mod.sigma = coef(summary(mod.1))[ , "Std. Error"]
#N.B. this can also be calculated from the model, and for each individual including individual effects
    # modmat  =  model.matrix(terms(mod.1), mdata)#extract model matrix
    # fixedfx = modmat %*% fixef(mod.1) #calulate fixed effects
    # pvar = diag(modmat %*% tcrossprod(vcov(mod.1), modmat)) #use variance-covariance matrix to calculate model variance
    # tvar = pvar+as.numeric(VarCorr(mod.1))
    # p_upper = fixedfx + 1.96*sqrt(tvar) 
    # p_lower = fixedfx - 1.96*sqrt(tvar)
#extract fixed-effects confidence interval (explained below)
#confidence interval
mod.fixef.ci = qnorm(c(0.05/2,1-0.05/2),
                    mean = mod.mu,
                    sd = mod.sigma
                    )
message('two-tailed 95% confidence interval for mean cross choice rate\n',
        round(plogis(mod.fixef.ci[1]), 3)*100,
        ' - ',
        round(plogis(mod.fixef.ci[2]), 3)*100,
        '%')

# Calculate p value -------------------------------------------------------
#critical p-value depends on the hypothesis
alpha = switch(EXPR = h1, #h1 should be a string with three possible values
               less = 1-0.05,#cumulative probability below mean should be 5%
               greater = 0.05,#cumulative probability _above_ mean should be 5%
               two_tailed = c(0,1)+c(1,-1)*0.05/2 #cumulative probability at either tail should be 5%
              )

# . Plot derived variables ------------------------------------------------
#Logistic regression assumes a normal distribution on the logit scale
#The fitted model has the probability density function
plot(x = seq(-5,5,0.01),
     y = dnorm(seq(-5,5,0.01),
           mean = mod.mu,
           sd = mod.sigma),
     type = 'l',
     xlab = 'log(odds)',
     ylab = 'probability density',
     main = 'fitted distribution for population'
     )
#For which the cumulative probability of the true mean
# being at a specific value, or lower, is
plot(x = plogis(seq(-5,5,0.01)),
     y = switch(EXPR = h1, #h1 should be a string with three possible values
                less = 1-pnorm( seq(-5,5,0.01),
                              mean = mod.mu,
                              sd = mod.sigma
                            ),#cumulative probability below mean should be 5%
                greater = pnorm( seq(-5,5,0.01),
                                   mean = mod.mu,
                                   sd = mod.sigma
                                   ),#cumulative probability _above_ mean should be 5%
                two_tailed = 2* sapply( seq(-5,5,0.01),
                                              function(x)
                                              {
                                        ifelse(test = x<mod.mu,
                                               yes = pnorm(q = x,
                                                          mean = mod.mu,
                                                          sd = mod.sigma 
                                                           ),
                                               no = 1-pnorm(q = x,
                                                            mean = mod.mu,
                                                            sd = mod.sigma 
                                                           )
                                              )
                                              }
                                        ), #cumulative probability at either tail should be 5%
               ),
     type = 'l',
     xlab = 'mean correct choice rate',#'log(odds)',
     ylab = switch(EXPR = h1, #h1 should be a string with three possible values
                   less = 'probability of true mean > value',#cumulative probability below mean should be 5%
                   greater = 'probability of true mean < value',#cumulative probability _above_ mean should be 5%
                   two_tailed = 'probability of true mean > or < value (maximum)' #cumulative probability at either tail should be 5%
     ),
     main = 'fitted cumulative probability'
     )
abline(h = switch(h1,
                  two_tailed = 2*min(alpha),
                  less = 1-alpha,
                  greater = alpha
                  ),
       col = 'red',
       lty = 2#dashed line
       )
#below the 5th percentile for this normal distribution (quantile function),
#the probability of a true mean at that value or lower is less than 0.05
abline(v = plogis(
                  qnorm(p = alpha,
                         mean = mod.mu,
                         sd = mod.sigma)
                  ),
       col = 'red',
       lty = 1#solid line
       )
#add null hypothesis choice rate
abline(v = h0,
       col = 'blue',
       lty = 1#solid line
       )
legend(x = 'topleft',
       legend = c('critical probability',
                  'critical choice rate',
                  'null hypothesis'
                  ),
       lty = c(2,1,1),
       col = c('red','red','blue'),
       cex = 0.35
      )
#Convert reference value to logit space
h0.logit = qlogis(h0)#quantile function for the logistic distribution performs logit tranform
#calculate z-score
z_score = (mod.mu-qlogis(h0))/mod.sigma
message('z-score = ',
        round(z_score,3)
        )
#Update, use z-score to calculate p-value
# p.h0 = ifelse(h1 == 'greater',#if alternative hypothesis that ture mean is greater than reference
#               pnorm(h0.logit,#cumulative probability up to reference
#                     mean = mod.mu,#distribution with model mean
#                     sd = mod.sigma#and model standard deviation
#                     ),
#               1-pnorm(h0.logit,#cumulative probability down to reference (1-p)
#                       mean = mod.mu,#distribution with model mean
#                       sd = mod.sigma#and model standard deviation
#                       )
#               )
p.h0 = switch(EXPR = h1,
              greater = pnorm(-z_score),
              less = pnorm(z_score),
              two_tailed = 2*pnorm(-abs(z_score))
              )
message(
        'Probability mean correct choices are NOT ',
        switch(EXPR = h1,
               greater = c(
                            h1,
                            ' than ',
                            h0
                            ),
               less = c(
                            h1,
                            ' than ',
                            h0
                            ),
               two_tailed = c(
                            'different from ',
                            h0
                            )
                ),
        '\n p = ',
        round(p.h0,5)
        )
#save result
write.table(x = cbind(h0 = h0,#save relevant info in a table
                      h1 = h1,
                      mean_choice = plogis(mod.mu),
                      ci_02.5 = plogis(mod.fixef.ci[1]),
                      ci_97.5 = plogis(mod.fixef.ci[2]),
                      z_score = z_score,
                      p = p.h0
                      ),
            file = file.path(dirname(path_file),#same place as original with
                             paste(sub(pattern='.csv',#similar name
                                       replacement = '',
                                       x = basename(path_file)
                                       ),
                                  'GLMM mean test.csv'#ending in mean test
                                  )
                             ),
            row.names = FALSE,#rows do not need names
            sep = csv_sep #Use same separator as original csv file
  
)

# Plot mean and confidence intervals --------------------------------------


#fixed effects model in probability space
#mean correct
barplot(height = plogis(mod.mu),
        width = 0.1,
        xlim = c(0,1),
        ylim = c(0,1),
        space = 0.5,
        main = paste0('p(mean',
                      ifelse(h1=='greater','≤','≥'),
                      'reference) = ',
                      round(p.h0,3)
                      )
        )
#reference correct choice rate
abline(h = h0,
       col = 'red'
       )
#Fixed effects 95% confidence interval for the mean
#N.B. This confidence interval is two-tailed, it extends to the 2.5th percentile
# whereas our test is one-tailed, its "confidence interval" extends to the 5th percentile
arrows(x0 = 0.1,
      x1 = 0.1,
      y0 = plogis(mod.fixef.ci[1]),
      y1 = plogis(mod.fixef.ci[2]),
      code = 3,
      angle = 90,
      length = 0.1
      )
#add mean back for appearances
points(x = 0.1,
       y = plogis(mod.mu),
       pch = 19#a small solid dot
       )
legend(x = 'topright',
       legend = c('reference',
                  'observed mean'
                  ),
       col = c('red', 'black'),
       lty = c(1, NA),
       pch = c(NA, 19),
       cex = 0.5
       )
