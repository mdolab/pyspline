subroutine cgns_test()

#ifdef USE_CGNS 

print *,'CGNS enabled in this build'
  
#else

print *,'CGNS not enabled in this build'

#endif

end subroutine cgns_test
