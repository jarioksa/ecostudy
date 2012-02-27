`poplogist`  <-
    function(t, r, N0, K)
{
    K/(1 + (K-N0)/N0 * exp(-r * t))
}
