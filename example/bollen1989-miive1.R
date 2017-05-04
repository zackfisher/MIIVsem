# Example 1

bollen1989a_model <- '

    Eta1 =~ y1 + y2  + y3  + y4  
    Eta2 =~ y5 + y6  + y7  + y8    
    Xi1  =~ x1 + x2 + x3 

    Eta1 ~ Xi1  
    Eta2 ~ Xi1 
    Eta2 ~ Eta1 

    y1   ~~ y5
    y2   ~~ y4
    y2   ~~ y6
    y3   ~~ y7
    y4   ~~ y8
    y6   ~~ y8 

  '

miive(model = bollen1989a_model, data = bollen1989a)
