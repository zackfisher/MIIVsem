#' Industrialization-Democracy Data
#'
#' A dataset from Bollen (1989) containing measures of political 
#' democracy and industrialization for 75 developing countries 
#' in 1960 and 1965. The variables are as follows:
#'
#' \itemize{
#'   \item y1. freedom of the press, 1960 
#'   \item y2. freedom of political opposition, 1960 
#'   \item y3. fairness of elections, 1960 
#'   \item y4. effectiveness of elected legislature, 1960
#'   \item y5. freedom of the press, 1965
#'   \item y6. freedom of political opposition, 1965
#'   \item y7. fairness of elections, 1965
#'   \item y8. effectiveness of elected legislature, 1965
#'   \item x1. natural log of GNP per capita, 1960
#'   \item x2. natural log of energy consumption per capita, 1960
#'   \item x3. arcsin of square root of percentage of labor force in industry, 1960
#' }
#'
#' @docType data
#' @keywords datasets
#' @name bollen1989a
#' @usage bollen1989a
#' @format A data frame with 75 rows and 9 variables
#' 
#' @examples 
#' 
#' \dontrun{
#'   model <- '
#'     Eta1 =~ y1 + y2  + y3  + y4  
#'     Eta2 =~ y5 + y6  + y7  + y8    
#'     Xi1  =~ x1 + x2 + x3 
#'     Eta1 ~ Xi1  
#'     Eta2 ~ Xi1 
#'     Eta2 ~ Eta1 
#'     y1   ~~ y5
#'     y2   ~~ y4
#'     y2   ~~ y6
#'     y3   ~~ y7
#'     y4   ~~ y8
#'     y6   ~~ y8 devtools::build_win()

#'   '
#' }
#' 
#' @references 
#'   Bollen, K. A. (1989). Structural equation models. 
#'   New York: Wiley-Interscience.
NULL

#' Union sentiment data
#'
#' A dataset from McDonald and Clelland (1984) reanalyzed by
#' Bollen (1989) containing data on union sentiment of
#' southern nonunion textile workers. 
#'
#' \itemize{
#'   \item deferenc. deference (submissiveness) to managers
#'   \item laboract. support for labor activism
#'   \item unionsen. sentiment towards unions
#'   \item yrsmill. log of years spent in textile mill
#'   \item age. centered age
#' }
#'
#' @docType data
#' @keywords datasets
#' @name bollen1989b
#' @usage bollen1989b
#' @format A data frame with 173 rows and 5 variables
#' 
#' @examples 
#' 
#' \dontrun{
#'  model <- '
#'    unionsen ~  deferenc + laboract + yrsmill
#'    deferenc ~  age
#'    laboract ~  age + deferenc
#'    yrsmill  ~~ age
#'  '
#' }
#' 
#' @references 
#'   Bollen, K. A. 1989. Structural Equations with Latent Variables. 
#'   New York: Wiley
#'   
#'   McDonald, A, J., & Clelland, D. A. (1984). Textile Workers 
#'   and Union Sentiment. Social Forces, 63(2), 502–521.
NULL

#' Attractiveness and academic ability
#'
#' This data comes from a study by Felson and Borhnstedt (1979)
#' of perceived attractiveness and academic ability in teenagers, 
#' sixth through ninth grade. The six variables are perception of 
#' academic ability (academic), perception of physical 
#' attractiveness (attract), grade point average (gpa),
#' height, weight, and a strangers' rating of attractiveness 
#' (rating). 
#'
#' \itemize{
#'   \item acad. 
#'   \item athl. 
#'   \item attract.
#'   \item gpa. 
#'   \item height. 
#'   \item weight. 
#'   \item rating. 
#' }
#'
#' @docType data
#' @keywords datasets
#' @name felson1979
#' @usage felson1979
#' @format A data frame with 209 rows and 7 variables
#' 
#' @examples 
#' 
#' \dontrun{
#'   model <-  '
#'		 acad    ~ gpa + attract
#'     attract ~ height + weight + rating + acad
#'   '
#' }
#' 
#' 
#' @references 
#'   Felson, R.B. & Bohrnstedt, G.W. (1979). "Are the good 
#'   beautiful or the beautiful good?" The relationship between 
#'   children's perceptions of ability and perceptions of 
#'   physical attractiveness. Social Psychology Quarterly, 
#'   42, 386–392.
NULL

#' Subjective class data
#'
#' The following data is from Bollen (1989) using data from
#' Kluegel et al. (1977). These data include  measures of actual 
#' income (inc) and occupational prestige (occ),  measures 
#' of respondents' subjective assessments of income (subinc), 
#' occupational prestige (subocc), and overall SES status (subgen).
#'
#' \itemize{
#'   \item occ. actual occupational prestige
#'   \item inc. actual income
#'   \item subocc. respondents' subjective assessments of prestige
#'   \item subinc. respondents' subjective assessments of income
#'   \item subgen. respondents' subjective assessments of SES status
#' }
#'
#' @docType data
#' @keywords datasets
#' @name bollen1989c
#' @usage bollen1989c
#' @format A data frame with 432 rows and 5 variables
#' 
#' @examples
#' 
#'\dontrun{
#'  model <-    '
#'		subinc  ~ inc + subocc
#'    subocc  ~ occ + subinc
#'    subgen  ~ subinc + subocc
#'		subinc ~~ subocc + subgen
#'    subocc ~~ subgen
#'		inc    ~~ occ
#'  '
#'}
#' 
#' @references 
#' Bollen, K. A. 1989. Structural Equations with Latent Variables. 
#' New York: Wiley
#'   
#' Kluegel, J. R., Singleton, R., & Starnes, C. E. (1977). 
#' Subjective Class Identification: A Multiple Indicator 
#' Approach. American Sociological Review, 42(4), 599–611. 
NULL

#' Perceived accessibility data
#'
#' Data come from a survey that was conducted in rural clusters of 
#' Tanzania in 1993. The goal was to collect information on the perceived 
#' accessibility of a specific family planning facility that serviced each 
#' cluster. Six informants were chosen: 3 female and 3 male. New informants 
#' were chosen for each cluster. Each informant was independently asked to 
#' rate the accessibility of the facility, and how easy it was to get to the
#' facility. More specifically the women informants were asked to rate how 
#' women of childbearing age perceived the accessibility and easiness and 
#' men were asked to rate how accessible and easy men perceived access to 
#' the clinic to be. Higher values indicate greater accessibility and ease 
#' of travel. The female informants' ratings are 1 to 3 and the male 
#' informants' ratings are 4 to 6. 
#'
#' \itemize{
#'   \item access1. 
#'   \item access2. 
#'   \item access3. 
#'   \item access4. 
#'   \item access5. 
#'   \item access6. 
#'   \item easy1. 
#'   \item easy2. 
#'   \item easy3. 
#'   \item easy4. 
#'   \item easy5. 
#'   \item easy6. 
#' }
#'
#' @docType data
#' @keywords datasets
#' @name bollen1996
#' @usage bollen1996
#' @format A data frame with 220 rows and 12 variables
#' 
#' @examples
#'
#'\dontrun{
#'  model <- ' 
#'     femaleAccess  =~ access1 + access2 + access3
#'     maleAccess    =~ access4 + access5 + access6
#'     femaleEasy    =~ easy1   + easy2   + easy3
#'     maleEasy      =~ easy4   + easy5   + easy6 
#'  '
#'}
#'
#' @references 
#' Bollen, K. A., Speizer, I. S., & Mroz, T. A. (1996). Family Planning 
#' Facilities in Rural Tanzania: His and Her Perceptions of Time and Distance. 
NULL


#' Reisenzein data
#'
#' This dataset comes from Reisenzein (1986). In this paper Reisenzein designed 
#' a randomized experiment to test Weiner's attribution-affect model of 
#' helping behavior. According to this theory, whether people help others is 
#' determined by their anger or sympathy. Anger and sympathy are affected by 
#' perceived controllability.  If the individuals have gotten into difficult 
#' situations as a result of their own controllable actions, then this 
#' negatively affects sympathy and positively affects anger of the potential 
#' helpers. The opposite holds if the situation seems beyond the individuals’ 
#' control. This data comes from an experiment that describes a person 
#' collapsing and lying on the floor of a subway.  Subjects were told that the 
#' person was either drunk (controllable situation) or ill 
#' (uncontrollable situation). This randomized story was intended to affect 
#' perceptions of controllability, and controllability in turn affected 
#' feelings of sympathy and anger.  Finally, sympathy should positively affect 
#' helping behavior while anger would negatively affect helping.   
#'
#' \itemize{
#'   \item Z1. Eliciting Situation
#'   \item Z2. How controllable, do you think, is the cause of the person's 
#'   present condition? (1 = not at all under personal control, 9 = completely 
#'   under personal control).
#'   \item Z3. How responsible, do you think, is that person for his present 
#'   condition? (1 = not at all responsible, 9 = very much responsible).
#'   \item Z4. I would think that it was the person's own fault that he is in 
#'   the present situation. (1 = no. not at all. 9 = yes, absolutely so).
#'   \item Z5. How much sympathy would you feel for that person? (1 = none at 
#'   all. 9 = very much).
#'   \item Z6. I would feel pity for this person. (1 = none at all, 9 = very 
#'   much).
#'   \item Z7. How much concern would you feel for this person? (1 = none al 
#'   all, 9 = very much).
#'   \item Z8. How angry would you feel at that person? (1 = not at all, 9 = 
#'   very much).
#'   \item Z9. How irritated would you feel by that person? (1 = not at all, 9 
#'   = very much).
#'   \item Z10. I would feel aggravated by that person. (1 = not at all, 9 = 
#'   very much so).
#'   \item Z11.	How likely is it that you would help that person? (1 = 
#'   definitely would not help. 9 = definitely would help).
#'   \item Z12.	How certain would you feel that you would help the person? 
#'   (1 = not at all certain. 9 = absolutely certain).
#'   \item Z13.	Which of the following actions would you most likely engage in? 
#'   1 = not help at all; 2 = try to alert other bystanders, but stay 
#'   uninvolved myself; 3 = try to inform the conductor or another official in 
#'   charge; 4 = go over and help the person to a seat; 5 = help in any way 
#'   that might be necessary, including if necessary first aid and/or 
#'   accompanying the person to a hospital.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name reisenzein1986
#' @usage reisenzein1986
#' @format A data frame with 138 rows and 13 variables
#' 
#'
#' @references 
#' Reisenzein, R. (1986). A Structural Equation Analysis of Weiner's 
#' Attribution-Affect Model of Helping Behavior. Journal of Personality 
#' and Social Psychology, 50(6), 1123–33.
NULL