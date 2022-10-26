"""
Function calculates probability of treatment assignment in a trial, targeting 
equal 1:1 allocation, according to Adjustable Biased Coin Design (ABCD).  

 author: BaldiAntognini A, Giovagnoli A. 
  title: A new 'biased coin design' for the sequential allocation of two 
         treatments. 
journal: Applied Statistics
   year: 2004
 volume: 53
 number: 4
  pages: 651-64
"""
function abcd(N1::Int64, N2::Int64; a::Float64)
    # probability function proposed by A. Baldi Antognini and A. Giovagnoli
    F(x, a) = abs(x)^a/(abs(x)^a + 1)

    # current treatment imbalance
    d = N1 - N2

    # calculating probability of tretament asssgnment, given imbalance
    prb = abs(d) <= 1 ? 0.5 : (d < -1 ? F(d, a) : 1-F(d, a))    

    return prb
end