# Example 3

bollen1989a_model_r <- '

    Eta1 =~ y1 + l2*y2  + l3*y3  + l4*y4  
    Eta2 =~ y5 + l2*y6  + l3*y7  + l4*y8    
    Xi1  =~ x1 + x2 + 0.5*x3 

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

miive(model = bollen1989a_model_r, data = bollen1989a)
