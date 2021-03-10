


data1 = yy_[:,0]
data2 = y_[:,0]
plt.hist2d(data1,data2,bins=100,norm=LogNorm(),cmap='jet')
plt.xlim([lstmin,lstmax])
plt.ylim([lstmin,lstmax])
plt.colorbar()
plt.title('SLSTR_N-AHI')
plt.plot([lstmin,lstmax],[lstmin,lstmax],'k--')
plt.grid()
plt.show()