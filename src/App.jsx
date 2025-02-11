import React, { useState, useRef } from "react";
import { BrowserRouter as Router, Route, Routes, Link } from "react-router-dom";
import MoleculeIndex from "./components/ChemicalEditor";
import './styles/main.css';

const App = () => {
  const [isLoggedIn, setIsLoggedIn] = useState(false);
  const [showModal, setShowModal] = useState(false);
  const [showDropdown, setShowDropdown] = useState(false);
  const [showRegisterModal, setShowRegisterModal] = useState(false);
  const [errorMessage, setErrorMessage] = useState("");  // Store error message
  const[uploadMessage, setUploadMessage] = useState("");

  const fileInputRef = useRef(null);
  

  // Handle file selection & upload
  const handleUploadClick = () => {
    fileInputRef.current.click();
  };
  const handleFileChange = async(event) => {
    const selectedFile = event.target.files[0];
    if (!selectedFile) {
      setUploadMessage("Please select a file first.");
      return;
    }
    const formData = new FormData();
    formData.append("file", selectedFile);
    try{
      const response = await fetch("http://localhost:3000/upload", {
        method: "POST",
        body: formData,
    });
    const result = await response.json();
    setUploadMessage(result.message);
  } catch (error) {
    setUploadMessage("An error occurred while uploading the file.try again later.");
  }
};



  // Use useRef to get input element values
  const emailRefLogin = useRef(null);
  const passwordRefLogin = useRef(null);
  const captchaRefLogin = useRef(null);
  const emailRefRegister = useRef(null);
  const passwordRefRegister = useRef(null);
  const repasswordRef = useRef(null);
  const captchaRefRegister = useRef(null);
  const avatarRef = useRef(null);

  const validateEmail = (email) => {
    const emailRegex = /^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$/;
    return emailRegex.test(email);
  };
  //8-16 characters, at least one one special character
  const validatePassword = (password) => {
    const passwordRegex = /^(?=.*[!@#$%^&*(),.?":{}|<>])[A-Za-z\d!@#$%^&*(),.?":{}|<>]{8,16}$/;
    return passwordRegex.test(password);
  };
  
  


  // Toggle login modal
  const toggleLoginModal = () => {
    setShowModal(!showModal);
    setShowDropdown(false);
    setShowRegisterModal(false);
    setErrorMessage("");  // Clear previous error messages
  };

  // Toggle register modal
  const toggleRegisterModal = () => {
    setShowRegisterModal(!showRegisterModal);
    setShowModal(false);
    setErrorMessage("");  // Clear previous error messages
  };

  // Toggle dropdown menu
  const toggleDropdown = () => {
    setShowDropdown(!showDropdown);
    setShowModal(false);
  };

  // Handle login
  const handleLogin = () => {
    const email = emailRefLogin.current.value;
    const password = passwordRefLogin.current.value;
    const captcha = captchaRefLogin.current.value;

    setErrorMessage("");  // Clear previous error messages

    // Validation checks
    if (!email) {
      setErrorMessage("Please enter a valid email address!");
      return;
    }
    if (!validateEmail(email)) {
      setErrorMessage("Please enter a valid email address!");
      return;
    }
    if (!password) {
      setErrorMessage("Please input the password!");
      return;
    }
    if (!captcha) {
      setErrorMessage("Please input the captcha!");
      return;
    }

    // If validation passes, login the user
    setIsLoggedIn(true);
    setShowModal(false);
  };

  // Handle registration
  const handleRegister = () => {
    const email = emailRefRegister.current.value;
    const password = passwordRefRegister.current.value;
    const repassword = repasswordRef.current.value;
    const captcha = captchaRefRegister.current.value;

    setErrorMessage("");  // Clear previous error messages

    // Validation checks
    if (!email) {
      setErrorMessage("Please enter a valid email address!");
      return;
    }
    if (!validateEmail(email)) {
      setErrorMessage("Please enter a valid email address!");
      return;
    }
    if (!password) {
      setErrorMessage("Please create the password!");
      return;
    }
    if (!validatePassword(password)) {
      setErrorMessage("Password must be 8-16 characters, at least one special character!");
      return;
    }
    if (!repassword) {
      setErrorMessage("Please confirm the password!");
      return;
    }
    if (repassword !== password) {
      setErrorMessage("Passwords don't match, please confirm!");
      return;
    }
    if (!captcha) {
      setErrorMessage("Please input the captcha!");
      return;
    }

    setIsLoggedIn(false);
    setShowRegisterModal(false);
  };

  // Handle logout
  const handleLogout = () => {
    setIsLoggedIn(false);
    setShowDropdown(false);
  };

  return (
    <Router>
      <Routes>
        <Route
          path="/"
          element={
            <div className="app-container">
              {/* Header */}
              <header className="app-header">
                <div className="user-info">
                  <button
                    className="avatar-button"
                    onClick={isLoggedIn ? toggleDropdown : toggleLoginModal} // Toggle login modal or user menu
                    ref={avatarRef} 
                  >
                    <img src="/assets/user.png" className="avatar" alt="User Avatar" />
                    {isLoggedIn ? <span className="avatar-name">User</span> : <span className="login-text">Login</span>}
                  </button>
                </div>
                <nav className="app-nav">
                  <a href="#">Structure</a>
                  <a href="#">Dictionary</a>
                  <a href="#">Cloud</a>
                  <a href="#">Message</a>
                </nav>
              </header>

              {/* Main page */}
              <main className="app-main">
                <h1 className="main-title">Create</h1>
                <section className="create-section">
                  <div className="action-buttons">
                    <Link to="/editor" className="action-button">
                      <div className="icon-placeholder">+</div>
                      <p>New Canvas</p>
                    </Link>
                    <button className="action-button" onClick={handleUploadClick}>
                      <div className="icon-placeholder">📂</div>
                      <p>Upload File</p>
                    </button>
                    <input
                      type = "file"
                      ref={fileInputRef}
                      style= {{display: "none"}}
                      onChange={handleFileChange}
                      />
                  </div>
                  {uploadMessage && <p className="upload-message">{uploadMessage}</p>}
                </section>

                {/* Recent Files */}
                <h2 className="main-title">Recent Files</h2>
                <section className="recent-files">
                  <div className="empty-state">
                    <img src="/path-to-empty-icon.png" alt="No Recent Files" />
                    <p>No recent files currently</p>
                  </div>
                </section>
              </main>

              {/* Login Modal */}
              {showModal && !isLoggedIn && !showRegisterModal && (
                <div className="modal-overlay" onClick={() => setShowModal(true)}>
                  <div className="modal" onClick={(e) => e.stopPropagation()}>
                    <button className="close-button" onClick={() => setShowModal(false)}>
                      &times;
                    </button>

                    <h2>SuperGlue</h2>
                    {errorMessage && (
                      <div className="error-message" 
                      style={{
                        color: "rgb(219, 33, 33)",
                        marginLeft: "0", 
                        marginBottom: "10px",  
                        maxWidth: "100%", 
                        wordWrap: "break-word", 
                        textAlign: "left", 
                      }}>
                        {errorMessage}
                      </div>
                    )}
                    <div className="input-group">
                      <span className="icon">📧</span>
                      <input type="text" ref={emailRefLogin} placeholder="Please enter your email" />
                    </div>
                    <div className="input-group">
                      <span className="icon">🔒</span>
                      <input type="password" ref={passwordRefLogin} placeholder="Please enter the password" />
                    </div>
                    <div className="input-group">
                      <span className="icon">🛡️</span>
                      <input type="text" ref={captchaRefLogin} placeholder="Input the captcha" />
                    </div>
                    <button className="login-button" onClick={handleLogin}>Log In</button>
                    <div className="new-user">
                      <button onClick={toggleRegisterModal} 
                      style={{
                        background: "none", 
                        border: "none", 
                        color: "#ffaa00", 
                        cursor: "pointer",
                        marginLeft: '199px',
                        marginTop: '25px'}}>
                        New User
                      </button>
                    </div>
                  </div>
                </div>
              )}

              {/* Register Modal */}
              {!showModal && !isLoggedIn && showRegisterModal && (
                <div className="modal-overlay" onClick={() => setShowRegisterModal(true)}>
                  <div className="modal-register" onClick={(e) => e.stopPropagation()}>
                  <button className="close-button2" onClick={() => setShowRegisterModal(false)}>
                      &times;
                    </button>
                    <h2>Register</h2>
                    {errorMessage && (
                      <div className="error-message" 
                      style={{
                        color: "rgb(219, 33, 33)",
                        marginLeft: "0", 
                        marginBottom: "10px",  
                        maxWidth: "100%", 
                        wordWrap: "break-word", 
                        textAlign: "left", 
                      }}>
                        {errorMessage}
                      </div>
                    )}
                    <div className="input-group">
                      <span className="icon">📧</span>
                      <input type="text" ref={emailRefRegister} placeholder="Please enter your email" />
                    </div>
                    <div className="input-group">
                      <span className="icon">🔒</span>
                      <input type="password" ref={passwordRefRegister} placeholder="Please create the password" />
                    </div>
                    <div className="input-group">
                      <span className="icon">🔒</span>
                      <input type="password" ref={repasswordRef} placeholder="Please confirm the password" />
                    </div>
                    <div className="input-group">
                      <span className="icon">🛡️</span>
                      <input type="text" ref={captchaRefRegister} placeholder="Input the captcha" />
                    </div>
                    <button className="register-button" onClick={handleRegister}>Register</button>
                    <div className="exsisting-user">
                      <button onClick={toggleLoginModal} 
                      style={{
                        background: "none", 
                        border: "none", 
                        color: "#ffaa00", 
                        cursor: "pointer"}}>
                        Already had an account?
                      </button>
                    </div>
                  </div>
                </div>
              )}

              {/* User Dropdown */}
              {isLoggedIn && showDropdown && (
                <div
                  className="dropdown-menu"
                  style={{
                    position: "absolute",
                    top: avatarRef.current.offsetTop + avatarRef.current.offsetHeight + 5 + "px",
                    left: avatarRef.current.offsetLeft + "px",
                  }}
                  onClick={(e) => e.stopPropagation()}
                >
                  <ul>
                    <li>Profile</li>
                    <li>Change Password</li>
                    <li>Feedback</li>
                    <li>Download App</li>
                    <li>About SuperGlue</li>
                    <li onClick={handleLogout}>Logout</li>
                  </ul>
                </div>
              )}
            </div>
          }
        />
        <Route path="/editor" element={<MoleculeIndex />} />
      </Routes>
    </Router>
  );
};

export default App;
