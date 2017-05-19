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