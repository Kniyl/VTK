// advanced.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// advanced dialog

class advanced : public CDialog
{
// Construction
public:
	advanced(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(advanced)
	enum { IDD = IDD_ADVANCED };
	CString	m_EXTRA_CFLAGS;
	CString	m_EXTRA_LINK_FLAGS;
	CString	m_WhereTcl;
	CString	m_WhereTk;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(advanced)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(advanced)
	afx_msg void OnWhereLibTCL();
	afx_msg void OnWhereLibtk();
	afx_msg BOOL Browse(CString& result, const char*);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
