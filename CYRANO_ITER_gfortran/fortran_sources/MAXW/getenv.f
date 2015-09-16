      integer function getenv(name, jobnam, i)

c     dummy function for compatibility with unicos version of Cyrano

      CHARACTER*20 NAME
      CHARACTER*7 JOBNAM
      integer I

      jobnam='Mac_Job'
      getenv = 0
      return
      end