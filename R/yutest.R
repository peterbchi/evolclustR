#' An evolClustR Function
#'
#' Function to return fractions in each cluster for any given radius
#' @param Cdata PDB file parsed to contain only central Carbon atoms
#' @param subs substited sites
#' @keywords cluster
#' @export
#' @examples
#' yu.test()

# comments from pg. 684
yu.test<-function(Cdata, subs, perm.reps){
  # For branch i, define Ri as the total number of amino acid
  # replacements on the branch and Ni as the average number of
  # replacements among the Ri 10-A balls corresponding to the branch.
  Ni <- mean(count.subs(Cdata, subs, radius=10))

  perm.null <- permute(Calphas, subs, reps=perm.reps, radius=10, contact=FALSE)

  # we use Nti
  # to be the average number of replacements per
  # ball on branch i for the tth permuted data set.
  Nti <- apply(apply(perm.null, 2, count.subs, Cdata=Calphas, radius=10), 2, mean)

  # With a total of T permuted matrices,
  T <- perm.reps

  # Equation (1):
  NiS <- sum(Nti) / T

  # Equation (2):
  sigmaiS2 <- (sum((Nti - NiS)^2)/(T-1))

  # Equation (3):
  Zi <- (Ni - NiS) / sqrt(sigmaiS2)

  # Rather than relying on any assumption of a specific form for the null distribution of Z
  # (e.g., a standard normal distribution), we approximate this distribution.
  # The approximation is straightforward to obtain by calculating
  # Z as in the procedure described above, except that Nti values
  # from permutations are substituted for Ni in Eq. 3. This generates T
  # simulated values of Z under the IS null hypothesis.
  Zti <- (Nti - NiS) / sqrt(sigmaiS2)

  # a p-value can be approximated by adding 1 to the number of simulated Z values
  # that exceed the observed value and then dividing by T+1.
  p.value <- (sum(Zti>Zi) + 1) / (T+1)

  return(p.value)

}

