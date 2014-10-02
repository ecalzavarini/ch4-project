<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <!-- <xsl:variable name="what">1</xsl:variable> -->
  <xsl:template name="treeview">
    <xsl:param name="depth"/> 
    <xsl:param name="content"/> 
    <xsl:param name="radix"/>
    <xsl:param name="radix_base"/>
    <xsl:param name="tot"/> 
    <xsl:param name="father"/> 

    <!--
	<xsl:for-each select="$content/variable">
	  ciao      <xsl:value-of select="@name"/>
	</xsl:for-each>
	-->

    <xsl:if test="not($content/core)">
      <xsl:variable name="id1"><xsl:value-of select="generate-id()"/></xsl:variable>
      <xsl:variable name="uniq_id">depth <xsl:value-of select="number($depth)"/> count a<xsl:value-of select="$id1"/></xsl:variable>
        [<xsl:value-of select="$uniq_id"/> father <xsl:value-of select="$father"/> ]</xsl:if>
    <!-- <xsl:value-of select="number($depth)"/> -->
    <xsl:for-each select="$content/core">
      <!-- print name on screen -->
      <!-- unique identifier of level-position -->
      
      <!-- <xsl:variable name="ia"><xsl:number/></xsl:variable> -->
      <xsl:variable name="ia"><xsl:value-of select="generate-id()"/></xsl:variable>
      <xsl:variable name="iddd"><xsl:value-of select="$father"/></xsl:variable>
      <xsl:variable name="id"><xsl:value-of select="$ia"/></xsl:variable>
      <xsl:variable name="uniq_id">depth <xsl:value-of select="number($depth)"/> count <xsl:value-of select="$id"/></xsl:variable>
      <xsl:variable name="radixc"><xsl:value-of select="$radix"/><xsl:value-of select="@name"/></xsl:variable>
      <xsl:if test="number($what)=4">
	[<xsl:value-of select="$uniq_id"/> father <xsl:value-of select="$father"/> ] <xsl:value-of select="translate($radixc,'_', '-')"/></xsl:if>
      <!-- <xsl:if test="number($id)=number($tot)">;</xsl:if> -->
      
      <xsl:for-each select="./variable">
	<!-- Variable --> <xsl:variable name="radixa"><xsl:value-of select="$radixc"/>_<xsl:value-of select="@name"/></xsl:variable>
	<xsl:if test="number($what)=4">
    	  <!-- Variable --> 
	  <xsl:variable name="aaa"><xsl:value-of select="@elements"/></xsl:variable>
	  - father <xsl:value-of select="$id"/> variable<xsl:text>  </xsl:text><xsl:value-of select="@type"/><xsl:text>  </xsl:text><xsl:value-of select="translate($radixa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/><xsl:text> </xsl:text> <xsl:value-of select="@value"/> <xsl:text> </xsl:text><xsl:value-of select="translate($aaa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>
	  <!-- <xsl:value-of select="$radixa"/> -->
	</xsl:if>
      </xsl:for-each>
      <xsl:choose>
	<xsl:when test="number($depth)=0">
	  <!-- <xsl:value-of select="@name"/> -->
	  <xsl:variable name="radixa"><xsl:value-of select="@name"/></xsl:variable>
	</xsl:when> 
	<xsl:otherwise>
	  <!-- <xsl:value-of select="$radix"/>_<xsl:value-of select="@name"/> -->
	  <xsl:variable name="radixa"><xsl:value-of select="$radix"/>_<xsl:value-of select="@name"/></xsl:variable>
	</xsl:otherwise>
      </xsl:choose>
      
      <xsl:variable name="radixa"><xsl:value-of select="$radix"/><xsl:value-of select="@name"/>_</xsl:variable>
      <xsl:variable name="radixb"><xsl:value-of select="$radix"/><xsl:value-of select="@name"/></xsl:variable>
      <!-- <xsl:value-of select="$radixa"/> -->
      <xsl:call-template name="treeview">
	<xsl:with-param name="content" select="."/>
	<xsl:with-param name="depth" select="number($depth)+1"/> 
	<xsl:with-param name="radix" select="$radixa"/>
	<xsl:with-param name="radix_base" select="$radixb"/>
	<xsl:with-param name="tot" select="$tot"/>
	<xsl:with-param name="father" select="$id"/>
      </xsl:call-template>
    </xsl:for-each>
    <!-- <xsl:value-of select="number(360) div number ($id)"/> -->
  </xsl:template>
  

  <xsl:template match="/">
    <xsl:variable name="elements"><xsl:value-of select="count(/lbe3d/core)"/></xsl:variable>
    \begin{tikzpicture}[mindmap,
    level 1 concept/.append style={level distance=200,sibling angle=<xsl:value-of select="number(360) div number ($elements)"/>},
    extra concept/.append style={color=blue!50,text=black}]
    \begin{scope}[mindmap, concept color=orange, text=black]
    \node [concept] {{\bf LBM3D}\\Federico Toschi}[clockwise from=-5]<xsl:call-template name="treeview">
      <xsl:with-param name="content" select="/lbe3d"/>
      <xsl:with-param name="depth" select="0"/>
      <xsl:with-param name="tot" select="$elements"/>
      <xsl:with-param name="father" select="1"/>
    </xsl:call-template>
    \end{scope}
    \end{tikzpicture}
  </xsl:template>
</xsl:stylesheet>
